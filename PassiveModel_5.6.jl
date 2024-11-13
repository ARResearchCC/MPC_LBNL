# Passive Model V5.6 Fred Fan 10/23/2024

############ Initialize Tools ############
begin
    import Pkg;
    # Initialize JuMP to allow mathematical programming models
    # Add Packages if you are running this for the first time
    
    # Add required packages
    # Pkg.add(["JuMP", "CSV", "DataFrames", "Dates", "XLSX", "LinearAlgebra", "StatsBase", "PyCall", "Statistics", "Clustering"])

    # Import necessary libraries
    using JuMP, CSV, DataFrames, Dates, XLSX, LinearAlgebra, StatsBase, PyCall, Statistics, Clustering
    include("Utility Functions.jl")
    include("Input_Parameters.jl")
end

############ Declare Parameters ############

begin
    # Define the columns and their bounds for normalization
    cols_for_clustering_SHGC = [:WS, :RH, :Ta, :Cloud, :GHI, :k_t]
    cols_for_clustering_RCC = [:WS, :RH, :Ta, :Cloud]
    bounds = Dict(
        :WS => (0.0, 10.0),  # Example bounds for Wind Speed # [m/s]
        :RH => (0.0, 100.0), # Example bounds for Relative Humidity # [%]
        :Ta => (0.0, 40.0), # Example bounds for Temperature in Celsius  # [°C]
        :GHI => (0.0, 1000.0), # Example bounds for Global Horizontal Irradiance [W/m^2]
        :Cloud => (0.0, 1.0), # Example bounds for Cloud cover fraction
        :k_t => (0.0, 1.0)  # Example bounds for clearness index
    )
end

function passive_model(file_path, df, Ti_constant, conversion_setting=1)
    begin
        #=
        Input raw weather data (Ta, Ti, WS, RH, GHI, cloudcover), tuned parameters
        1. Convert all weather data to SI
        2. Assume T_indoor is maintained at Ti_constant [°C], calculate Q_cd, Q_cv, Q_rh, and Q_rc [W]
        3. Report back Q_total and Q_cd, Q_cv, Q_rh, and Q_rc in [W]
        =#
    end
    
    begin
        # Load the Excel file
        # file_path = "C:\\Users\\Fred\\Desktop\\Microgrid_Project\\Test_Results\\Calibration_stepsize=5._V5.6._2024-10-22\\Calibration_Model5.6)_2024-10-22.2024-10-22_17-41-20.xlsx"
    
        # Create a dictionary to hold DataFrames for each sheet
        dfs = Dict{String, DataFrame}()
    
        # Open the Excel file
        XLSX.openxlsx(file_path) do workbook
            # Get the names of all sheets
            sheet_names = XLSX.sheetnames(workbook)
            
            # Loop through each sheet name and read it into a DataFrame
            for sheet_name in sheet_names
                sheet = XLSX.gettable(workbook[sheet_name])  # Read each sheet as a table
                dfs[sheet_name] = DataFrame(sheet)           # Convert to DataFrame
            end
        end
        
        SHGCs = dfs["Tuned Parameters"][:, "_SHGC_"]
        RCCs_day = dfs["Tuned Parameters"][:, "_RCC_day_"]
        RCCs_night = dfs["Tuned Parameters"][:, "_RCC_night_"]
        TC = dfs["Other Parameters"][1, "_Thermal_Capacitance_"]
        
        day_centroids_SHGC_df = dfs["SHGC Day Centroids"]
        day_centroids_RCC_df = dfs["RCC Day Centroids"]
        night_centroids_RCC_df = dfs["RCC Night Centroids"]
        
        day_centroids_SHGC = permutedims(Matrix{Float64}(day_centroids_SHGC_df))
        day_centroids_RCC = permutedims(Matrix{Float64}(day_centroids_RCC_df))
        night_centroids_RCC = permutedims(Matrix{Float64}(night_centroids_RCC_df))
    end

    NumTime = size(df)[1]

    # Create space to store heat transfers
    NetHeatTransfers = [] # [W]
    
    # Create space to store various sources of heat transfers.
    Q_cond = [] # [W]
    Q_conv = [] # [W]
    Q_radcool = [] # [W]
    Q_radheat = [] # [W]

    TemperatureAmbient = df[:, :Ta] # [°C]
    # TemperatureIndoor = df[:, :Ti] # [°C]
    WindSpeed = df[:, :WS] # [m/s]
    RelativeHumidity = df[:, :RH] # [%]
    GHI = df[:, :GHI] # [W/m^2]

    local_time = df[:, :DateTime]
    k_t = df[!, :k_t]
    solar_time = df[!, :SolarTime]
    CloudCover = df[!, :Cloud]
    day_status = df[!, :DayStatus]

    # assign cluster based on the distance closest to the centroid

    # Timeseries simulation
    for t = 1:NumTime
        # Conduction
        q_cd = Q_conduction(TemperatureAmbient[t], Ti_constant)/1000 # [kW]
        push!(Q_cond, q_cd)
        # Convection
        q_cv = Q_convection(WindSpeed[t], TemperatureAmbient[t], Ti_constant)/1000 # [kW]
        push!(Q_conv, q_cv)
        
        # Radiative cooling
        q_rc = Q_radiativeCooling(Ti_constant, TemperatureAmbient[t], RelativeHumidity[t], GHI[t], solar_time[t], day_status[t])/1000 # [kW]

        # Radiative heat gain to the structure
        q_rh = calculate_solarheatgain(solar_time[t], GHI[t], k_t[t], local_time[t])/1000 # [kW]
        
        new_row = DataFrame(WS = WindSpeed[t], RH = RelativeHumidity[t], Ta = TemperatureAmbient[t], GHI = GHI[t], Cloud = CloudCover[t], k_t = k_t[t])
        if day_status[t] == 0.0
            normalized_new_row_rcc = normalize_row(new_row, cols_for_clustering_RCC, bounds)
            # Find the nearest centroid for night time data
            nearest_rcc_cluster = find_nearest_centroid(normalized_new_row_rcc, night_centroids_RCC, cols_for_clustering_RCC)
            push!(Q_radcool, q_rc * RCCs_night[nearest_rcc_cluster])
            push!(Q_radheat, 0)
            Q_net = q_cd + q_cv - q_rc * RCCs_night[nearest_rcc_cluster] # [kW]
            push!(NetHeatTransfers, Q_net) 
        else
            normalized_new_row_shgc = normalize_row(new_row, cols_for_clustering_SHGC, bounds)
            normalized_new_row_rcc = normalize_row(new_row, cols_for_clustering_RCC, bounds)
            # Find the nearest centroid for day time data
            nearest_shgc_cluster = find_nearest_centroid(normalized_new_row_shgc, day_centroids_SHGC, cols_for_clustering_SHGC)
            nearest_rcc_cluster = find_nearest_centroid(normalized_new_row_rcc, day_centroids_RCC, cols_for_clustering_RCC)
            push!(Q_radcool, q_rc * RCCs_day[nearest_rcc_cluster])
            push!(Q_radheat, q_rh * SHGCs[nearest_shgc_cluster])

            Q_net = q_cd + q_cv + q_rh * SHGCs[nearest_shgc_cluster] - q_rc * RCCs_day[nearest_rcc_cluster] # [kW]
            push!(NetHeatTransfers, Q_net)
        end
    end    

    return NetHeatTransfers, Q_cond, Q_conv, Q_radheat, Q_radcool, TC
end
