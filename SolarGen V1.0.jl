# Solar Generation V1.0
# Complete Date: 10/28/2024
# ACE Microgrid Project --- CEE Atmosphere and Energy --- Stanford University

# This code is an intellectual property of Yuanbei "Fred" Fan.
# Dear user, this is Fred. If there are any questions with regards to the code, please do not feel hesistant to reach out to me at yf1098@stanford.edu.

############ Initialize Tools ############
begin
    using Pkg

    # Please unmute the next two lines if you are running this code for the first time and don't have all the packages installed.
    # Pkg.add("JuMP", "CSV", "DataFrames", "Dates", "Statistics", "PyCall")
    
    using JuMP
    using CSV
    using DataFrames
    using Dates
    using Statistics
    using PyCall

    # Code setup
    ENV["PYTHON"] = "C:\\Users\\Fred\\anaconda3\\python.exe"
    include("Input_Parameters.jl")

    # Import necessary Python libraries
    pd = pyimport("pandas")
    pvlib = pyimport("pvlib")
end

function Generate_PV(weather)
    
    # Convert Julia DataFrame to Pandas DataFrame
    columns = names(weather)
    data = [getproperty(weather, col) for col in columns]
    pandas_df = pd.DataFrame(Dict(zip(columns, data)))

    # Ensure the datetime column is in the correct format and set as the index
    pandas_df["DateTime"] = pd.to_datetime(pandas_df["DateTime"])
    pandas_df = pandas_df.set_index("DateTime")
    
    # Calculate PV module cell temperature with PVLib [°C]
    begin
        # Prepare data input for pvlib cell temperature calculation
        ambient_temperature = pandas_df["Ta"] # [°C]
        wind_speed = pandas_df["WS"] # [m/s]
        total_irradiance = pandas_df["GHI"] # [W/m^2]

        # Calculate PV cell temperature using the Fuentes Model
        temp_fuentes = pvlib.temperature.fuentes(
            total_irradiance,
            ambient_temperature,
            wind_speed,
            noct_installed,
            module_height=module_height,
            wind_height=wind_height,
            emissivity=module_emissivity,
            absorption=module_absorption,
            surface_tilt=module_surface_tilt,
            module_width=module_width,
            module_length=module_length
        ) # [°C]

        # Store PV cell temperature into weather file
        temp_fuentes_array_manual = [temp_fuentes[i] for i in 1:length(temp_fuentes)]
        weather.celltemp = temp_fuentes_array_manual
    end
    
    Num_Time = size(weather)[1]

    # Calculate PV capacity factor [kW DC Output/kW Capacity] 
    PV_CF = zeros(Num_Time)
    for i = 1:Num_Time
        pv_o = calculate_PV(weather[i, :DateTime], weather[i, :DNI], weather[i, :DHI], weather[i, :GHI], 180-Berg_tilt, module_surface_tilt, weather[i, :celltemp])
        PV_CF[i] = pv_o
    end

    return PV_CF # [kW/kW capacity]
end

# Define a function to calculate PV generation (capacity factor) for each time instance 
function calculate_PV(local_time, DNI, DHI, GHI, surface_azimuth, tilt_angle, T_cell)
    # Convert local time to Local Solar Time
    LT = hour(local_time) + minute(local_time)/60 + second(local_time)/3600

    # Day of the year (N)
    N = dayofyear(local_time)

    # Calculate the Equation of Time (Eqt)
    if 1 ≤ N ≤ 106
        Eqt = -14.2 * sin(π * (N+7)/111)
    elseif 107 ≤ N ≤ 166
        Eqt = 4 * sin(π * (N-106)/59)
    elseif 167 ≤ N ≤ 246
        Eqt = 6.5 * sin(π * (N-166)/80)
    elseif 247 ≤ N ≤ 366 #leap year
        Eqt = 16.4 * sin(π * (N-247)/113)  
    end

    # Calculate the Solar Time (T_solar)
    T_solar = LT + Eqt/60 + (4*(standard_meridian_longitude - longitude))/60 # [DEGREES]

    # Solar Declination (δ)
    δ = 23.45 * sin(deg2rad((360/365) * (284 + N))) # [DEGREES]

    # Solar Time Angle (ω)
    ω = 15 * (T_solar - 12) # [DEGREES]

    # Latitude (φ)
    φ = deg2rad(latitude) # [RADIANS]

    # Solar Elevation Angle (α)
    α = asin(sin(φ) * sin(deg2rad(δ)) + cos(φ) * cos(deg2rad(δ)) * cos(deg2rad(ω))) # [RADIANS]
    
    # Solar Azimuth Angle (Ψ)
    if ω < 0
        Ψ = acos((sin(deg2rad(δ)) - sin(α) * sin(φ)) / (cos(α) * cos(φ))) # [RADIANS]
    else
        Ψ = 2 * π - acos((sin(deg2rad(δ)) - sin(α) * sin(φ)) / (cos(α) * cos(φ))) # [RADIANS]
    end

    # Surface tilt (β)
    β = deg2rad(tilt_angle) # [RADIANS]
    # Surface azimuth (γ)
    γ = deg2rad(surface_azimuth) # [RADIANS]
    # Angle of Incidence (θ)
    θ = acos(sin(α) * cos(β) + cos(α) * sin(β) * cos(γ - Ψ)) # [RADIANS]

    # Beam Irradiance (I_b)
    # no direct beam irradiance if the angle of incidence is negative
    if cos(θ) < 0
        I_b = 0
    else
        I_b = DNI * cos(θ) # [W/m^2]
    end
    
    # Diffuse irradiance (I_d)
    I_d = max(DHI * (1+cos(β))/2 + 0.5*GHI*(0.012*(π/2 - α)-0.04)*(1-cos(β)),0) # [W/m^2]
    
    # Reflected irradiance (I_r)
    I_r = max(GHI * albedo * (1 - cos(β)) / 2, 0) # [W/m^2]
    
    # Total irradiance (I_poa)
    I_poa = I_d + I_r + I_b # [W/m^2]

    # Calculate transmittence 
    begin
        # Calculate angle of refraction into AR coating (θ_2) 
        θ_2 = asin((n_air / n_AR) * sin(θ)) # [RADIANS]
        
        # Calculate transmittance through AR coating (τ_AR) 
        τ_AR = 1 - 0.5 * ((sin(θ_2 - θ)^2) / (sin(θ_2 + θ)^2)) + ((tan(θ_2 - θ)^2) / ((tan(θ_2 + θ)^2))) # [1]

        # Calculate angle of refraction into glass cover (θ_3) 
        θ_3 = asin((n_AR / n_glass) * sin(θ_2)) # [RADIANS]

        # Calculate transmittance through glass cover (τ_glass) 
        τ_glass = 1 - 0.5 * ((sin(θ_3 - θ_2)^2) / (sin(θ_3 + θ_2)^2)) + ((tan(θ_3 - θ_2)^2) / ((tan(θ_3 + θ_2)^2))) # [1]
        
        # Calculate effective transmittance through AR coated module (τ_cover)
        τ_cover = τ_AR * τ_glass # [1]
    end
    
    # Calculate transmitted POA irradiance [I_tr]
    I_tr = I_poa * τ_cover # [W/m^2]

    # Calculate PV DC power output [P_dc]
    P_dc = (I_tr/1000) * P_dc0 * (1 + Γ_t * (T_cell - T_ref)) # [kW/kW]

    # Calculate PV DC power output after system loss [P_PV]
    P_PV = P_dc * η_PV # [kW/kW]

    return P_PV 
end 