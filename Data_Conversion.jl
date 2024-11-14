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

# Convert the pandas DataFrame to a Julia DataFrame
function pandas_to_julia(py_df)
    columns = py_df.columns.tolist()
    data_dict = Dict{String, Any}()
    for col in columns
        data_dict[col] = py_df[col].values
    end
    df = DataFrame(data_dict)
    select!(df, :Time, Not(:Time))
    return df
end

function resample_dataframe(df::DataFrame, current_step::Int, new_step::Int)
    # Check if the new step size is a multiple of the current step size
    if new_step % current_step != 0
        error("New step size must be a multiple of the current step size (5 minutes).")
    end

    # Calculate the ratio of the new step size to the current step size
    step_ratio = new_step ÷ current_step

    # Create a new DataFrame with the resampled data
    new_df = DataFrame()

    # Assuming the first column is `DateTime` and others follow
    new_df.DateTime = df.DateTime[1:step_ratio:end]

    # Copy and aggregate or sample other columns as needed
    for col_name in names(df)[2:end]
        # Placeholder: here you could perform averaging, summing, or other aggregations
        # This example simply takes the first value from each block of `step_ratio` rows
        new_df[!, col_name] = [df[!, col_name][i] for i in 1:step_ratio:nrow(df)]
    end

    return new_df
end

function prepare_df(df_raw)

    df_raw[!, :Ti] = (df_raw[!, :Ti] .- 32) .* 5/9 # [°C]
    df_raw[!, :Ta] = (df_raw[!, :Ta] .- 32) .* 5/9 # [°C]
    df_raw[!, :WS] = df_raw[!, :WS]  .* 0.44704 # [m/s]

    # Adding Columns
    rename!(df_raw, :Time => :DateTime)
    df_raw = create_ghi_ratio(df_raw)
    df_raw = add_solar_time(df_raw)
    df_raw = add_day_status(df_raw)

    return df_raw
end

function prepare_NSRDB(df_raw)

    rename!(df_raw, :Temperature => :Ta)
    rename!(df_raw, "Wind Speed" => :WS)
    rename!(df_raw, "Relative Humidity" => :RH)

    # Create the DateTime column
    df_raw.DateTime = DateTime.(df_raw.Year, df_raw.Month, df_raw.Day, df_raw.Hour, df_raw.Minute)
    df_raw = create_ghi_ratio(df_raw)
    df_raw = add_solar_time(df_raw)
    df_raw = add_day_status(df_raw)
    df_raw = add_cloudcover(df_raw)
    return df_raw
end

function prepare_epw(df_raw)

    rename!(df_raw, :Temperature => :Ta)
    rename!(df_raw, "Wind Speed" => :WS)
    rename!(df_raw, "Relative Humidity" => :RH)
    rename!(df_raw, "Dew Point" => :Td)

    # Create the DateTime column
    # df_raw.DateTime = DateTime.(df_raw.Year, df_raw.Month, df_raw.Day, df_raw.Hour, df_raw.Minute)
    rename!(df_raw, :timestamp => :DateTime)
    df_raw = create_ghi_ratio(df_raw)
    df_raw = add_solar_time(df_raw)
    df_raw = add_day_status(df_raw)
    df_raw = add_cloudcover_epw(df_raw)
    return df_raw
end

function convert_weather(df, input_format)
    if input_format == "site"
        df = prepare_df(df)
    elseif input_format == "NSRDB"
        df = prepare_NSRDB(df)
    elseif input_format == "epw"
        df = prepare_epw(df)
    end
    return df
end

function resample_mpc(df::DataFrame, start_index::Int, steps::Vector{Int}, stepsize::Int)
    # Initialize an empty DataFrame with the same columns as `df`
    selected_df = DataFrame([Symbol(col) => Any[] for col in names(df)])

    current_index = start_index
    for step in steps
        next_index = current_index + Int(step / stepsize)
        if next_index <= nrow(df)
            # Append the selected row to `selected_df`
            push!(selected_df, df[next_index, :])
            current_index = next_index
        else
            break
        end
    end

    return selected_df
end

function aggregate_energy(df::DataFrame, start_index::Int, steps::Vector{Int}, stepsize)
    selected_df = DataFrame()

    current_index = start_index

    QPassive = zeros(size(steps)[1])
    Lighting = zeros(size(steps)[1])
    Plugs = zeros(size(steps)[1])
    Occupancy = zeros(size(steps)[1])
    PV = zeros(size(steps)[1])

    for i = 1:length(steps)
        next_index = current_index + Int(steps[i] / stepsize)
        if next_index <= nrow(df)

            QPassive[i] = sum(df[current_index:(next_index-1), :"HVAC"]) * (stepsize/60) # [kWh]
            Lighting[i] = sum(df[current_index:(next_index-1), :"Lighting"]) * (stepsize/60) # [kWh]
            Plugs[i] = sum(df[current_index:(next_index-1), :"Plugs"]) * (stepsize/60) # [kWh]
            Occupancy[i] = sum(df[current_index:(next_index-1), :"Occupancy"]) * (stepsize/60) # [kWh]
            PV[i] = sum(df[current_index:(next_index-1), :"PV"]) * (stepsize/60) # [kWh/kW capacity]
            
            current_index = next_index
        else
            break
        end
    end

    selected_df.PV = PV
    selected_df.Lighting = Lighting
    selected_df.Plugs = Plugs
    selected_df.Occupancy = Occupancy
    selected_df.QPassive = QPassive

    return selected_df
end

function combine_dfs(df, Q_total, pv_cf, e_load)
    df.HVAC = Q_total
    df.PV = pv_cf
    df_new = hcat(df, e_load)
    return df_new
end

function load_generation(df, format)
    new_df = DataFrame()
    new_df.DateTime = df.DateTime
    new_df.PV = df.PV # [kW/kW capacity]
    # combine Q and E? what about the COP of HVAC system? It is not constant. (format = 2 means we assume constant COP H and C)
    
    if format == 1
        new_df.Q_load = df.HVAC .+ df.Occupancy .+ df.Lighting .+ df.Plugs # [kW]
        new_df.E_load = df.Lighting .+ df.Plugs # [kW]
    else
        E_load = zeros(size(new_df)[1])
        for i = 1:size(new_df)[1]
            # if the net Q is larger than 0, then we need a cooling load. If the net Q is smaller than 0, then we need a heating load.
            if (df.HVAC[i] + df.Occupancy[i] + df.Lighting[i] + df.Plugs[i]) >= 0
                E_load[i] = (df.HVAC[i] + df.Occupancy[i] + df.Lighting[i] + df.Plugs[i])/COP_C + df.Lighting[i] .+ df.Plugs[i] # [kW]
            else
                E_load[i] = (df.HVAC[i] + df.Occupancy[i] + df.Lighting[i] + df.Plugs[i])/COP_H + df.Lighting[i] .+ df.Plugs[i] # [kW]
            end
        end
        new_df.E_load = E_load # [kW] 
    end
    return new_df
end

function p2e(M_States)
    
    # Length of the M_States dataframe
    msize = size(M_States)[1];
    
    time_intervals = zeros(msize-1);

    for i = 1:msize-1
        time_intervals[i] = (M_States[i+1, :"Time"] - M_States[i, :"Time"])/3600 # [hr]
    end

    E_Pump_1_0 = sum(((M_States[i+1, :"Pump 1 Electric Power"] + M_States[i, :"Pump 1 Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]

    E_Pump_2_0 = sum(((M_States[i+1, :"Pump 2 Electric Power"] + M_States[i, :"Pump 2 Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]

    E_HP_0 = sum(((M_States[i+1, :"HP Electric Power"] + M_States[i, :"HP Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]

    Thermal_Power_Delivered_0 = sum(((M_States[i+1, :"HP Useful Thermal Power"] + M_States[i, :"HP Useful Thermal Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]

    E_Fan_0 = sum(((M_States[i+1, :"Fan Coil Fan Electric Power"] + M_States[i, :"Fan Coil Fan Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]

    E_PCM_C_0 = sum(((M_States[i+1, :"Heat Transfer Cold PCM"] + M_States[i, :"Heat Transfer Cold PCM"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]

    E_PCM_H_0 = sum(((M_States[i+1, :"Heat Transfer Hot PCM"] + M_States[i, :"Heat Transfer Hot PCM"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]

    E_FCU_0 = sum(((M_States[i+1, :"FCU Heat Delivered"] + M_States[i, :"FCU Heat Delivered"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
    
    E = [E_Fan_0, E_FCU_0, E_HP_0, E_PCM_C_0, E_PCM_H_0, E_Pump_1_0, E_Pump_2_0, Thermal_Power_Delivered_0]
            
    return E
end 

function find_average(M_States, col_name)
    
    # Length of the M_States dataframe
    msize = size(M_States)[1];
    
    time_intervals = zeros(msize-1);

    for i = 1:msize-1
        time_intervals[i] = (M_States[i+1, :"Time"] - M_States[i, :"Time"]) # [s]
    end

    time_length = M_States[msize, :"Time"] - M_States[1, :"Time"] # [s]
    
    average_value = (sum(((M_States[i+1, "$col_name"] + M_States[i, "$col_name"])/2)*time_intervals[i] for i = 1:(msize-1)))/time_length        
    return average_value
end 