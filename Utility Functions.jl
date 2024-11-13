# Utility Functions matched with Passive Model V5.6 Fred Fan 10/23/2024
############ Initialize Tools ############
begin
    import Pkg;
    # Initialize JuMP to allow mathematical programming models
    # Add Packages if you are running this for the first time
    
    # Add required packages
    # Pkg.add(["JuMP", "CSV", "DataFrames", "Dates", "XLSX", "LinearAlgebra", "StatsBase", "PyCall", "Statistics", "Clustering"])

    # Import necessary libraries
    using JuMP, CSV, DataFrames, Dates, XLSX, LinearAlgebra, StatsBase, PyCall, Statistics, Clustering
    include("Input_Parameters.jl")
end

function add_cloudcover(df)
    cloudcovers = zeros(size(df)[1])
    for i = 1:size(df)[1]
        cloudcovers[i] = cloud_type_to_percent(df[i, "Cloud Type"])
    end
    df.Cloud = cloudcovers
    return df
end    

function cloud_type_to_percent(cloud_type::Int)
    # Ensure cloud_type is within valid range
    if cloud_type < 0 || cloud_type > 9
        error("Invalid cloud type. Must be an integer between 0 and 9.")
    end

    # Map NSRDB cloud types to percentage cloud cover
    cloud_cover = Dict(
        0 => 0,   # Clear sky
        1 => 0.2,  # Cumulus - scattered clouds
        2 => 0.4,  # Altocumulus - partial cloud cover
        3 => 0.1,  # Cirrus - minimal obstruction
        4 => 0.25,  # Cirrostratus - thin clouds, partial coverage
        5 => 0.7,  # Stratus - thick cloud cover
        6 => 0.9,  # Nimbostratus - very thick, overcast, precipitation likely
        7 => 0.6,  # Stratocumulus - broken but significant cloud cover
        8 => 1, # Cumulonimbus - full cloud cover with storm potential
        9 => 0.5   # Indeterminate - assumed average cloud cover
    )

    return cloud_cover[cloud_type]
end

function add_cloudcover_epw(df) 
    cloudcovers = zeros(size(df)[1])
    for i = 1:size(df)[1]
        cloudcovers[i] = estimate_cloudcover(df[i, "Ta"], df[i, "Td"], df[i, "RH"])
    end
    df.Cloud = cloudcovers
    return df
end 

function create_ghi_ratio(df)

    start_date = df[1, :DateTime] + Hour(8)
    end_date = df[end, :DateTime] + Hour(8)
    stepsize = df[2, :DateTime] - df[1, :DateTime]
    ghi = calculate_ghi_julia(latitude, longitude, start_date, end_date, timezone, stepsize)
    df[!, :GHI] = max.(0, df[!, :GHI])
    df[!, :GHI_time] = ghi[:, 1]
    df[!, :GHI_perfect] = ghi[:, 2]
    df[!, :k_t] = zeros(size(df)[1])
    for i = 1:size(df)[1]
        if df[i, :GHI_perfect] > 10
            df[i, :k_t] = min(df[i, :GHI]/df[i, :GHI_perfect], 1)
        end
    end
    return df
end    

function correct_daylightsaving(df) 
    for i = 1:size(df)[1]
        if dls_start <= df[i, :DateTime] <= dls_end
            df[i, :DateTime] = df[i, :DateTime] - Hour(1)
        end
    end
    return df
end

function calculate_ghi_julia(lat, lon, start_date, end_date, timezone, stepsize)
    # Import necessary Python modules using PyCall
    pd = pyimport("pandas")
    pvlib = pyimport("pvlib")
    # Create a date range in Julia
    julia_dates = start_date:stepsize:end_date
    # Convert Julia dates to string format and then to Python list
    date_strings = [Dates.format(dt, "yyyy-mm-dd HH:MM:SS") for dt in julia_dates]
    python_dates = PyObject(date_strings)
    
    # Convert the list of date strings to a pandas DateTimeIndex with the specified timezone
    times = pd.to_datetime(python_dates).tz_localize("UTC").tz_convert(timezone)
    # times = pd.to_datetime(python_dates).tz_localize(timezone)
    # Create a location object in Python's pvlib
    location = pvlib.location.Location(lat, lon)
    
    # Calculate the solar position and clear sky data
    solar_position = location.get_solarposition(times)
    clearsky = location.get_clearsky(times)

    # Filter out nighttime values based on solar elevation using Python's True
    clearsky["ghi"].where(solar_position["elevation"] > 0, 0, inplace=PyObject(true))

    # Format the index in the desired "24:00:00" format for midnight
    formatted_index = [dt.strftime("%Y-%m-%d %H:%M:%S") for dt in clearsky.index.to_pydatetime()]

    # Create a DataFrame
    df = DataFrame(Time = formatted_index, GHI = clearsky["ghi"].to_numpy())
    return df
end

function add_day_status(df)
    day_status_list = zeros(size(df)[1])
    for i = 1:size(df)[1]
        solar_time = df[i, :SolarTime]
        local_time = df[i, :DateTime]
        day_status_list[i] = IsDay(solar_time, local_time)
    end
    df[!, :DayStatus] = day_status_list
    return df
end

function calculate_dew_point(T, RH)

    # Calculate the saturation vapor pressure (es) at the temperature T
    es = 6.112 * exp((17.67 * T) / (T + 243.5))
    
    # Calculate the actual vapor pressure (e)
    e = (RH / 100) * es
    
    # Calculate the dew point temperature (Td)
    Td = (243.5 * log(e / 6.112)) / (17.67 - log(e / 6.112)) # [°C]
    
    return Td
end

function Q_convection(WS, Ta, Ti)
    # calculate Q_convection using wind speed, ambient temperature, and indoor temperature
    # Method from 2001 ASHRAE HANDBOOK Page 2621
    
    # Convert temperatures from [°C] to [°K]
    T_indoor = Ti + 273.15 # [°K]
    ΔT = abs(T_indoor - T_d) # [°K]

    # Calculate infiltration
    # Units conversion has been made
    Infiltration = Al*sqrt(Cs*ΔT + Cw*WS^2) # [L/s]
    # Convert infiltration from [L/s] to [ft^3/min] (1[L/s] * 0.0353147[ft^3/L] * 60[s/min] = 2.11888[ft^3/min])
    # Convert infiltration from [L/s] to [m^3/s] 1[L/s]=0.001[m^3/s]
    #CFM = Infiltration * 2.11888 # [ft^3/min]  
    CFM = Infiltration * 0.001 # [m^3/s]
    # Calculate the heat exchange during infiltration (1.08 = 0.24[BTU/(lb*°F)] * 0.075[lb/ft^3] * 60[min/hr] (specific heat capacity of air at STP * air density at STP * min/hr conversion factor)
    
    # Q_inf_sen = (Ta - Ti) * 1.08 * CFM # [BTU/hr]
    
    # specific heat capacity of air = 1005 J/kg-K
    # air density at ATP = 1.25 kg/m^3
    Q_inf_sen = 1005 * 1.25 * CFM * (Ta - Ti) # [W]

    return Q_inf_sen
end

function Q_conduction(Ta, Ti)
    Q_cond = UA * (Ta - Ti) # [W] 
    return Q_cond
end

function add_solar_time(df)
    solartimes = zeros(size(df)[1])
    for i = 1:size(df)[1]
        local_time = df[i, :DateTime]
        # Convert local time to Local Solar Time
        LT = hour(local_time) + minute(local_time)/60 + second(local_time)/3600

        # Day of the year
        N = dayofyear(local_time)

        # Calculate the Equation of Time (Eqt)
        if 1 ≤ N ≤ 106
            Eqt = -14.2 * sin(π * (N+7)/111)
        elseif 107 ≤ N ≤ 166
            Eqt = 4 * sin(π * (N-106)/59)
        elseif 167 ≤ N ≤ 246
            Eqt = 6.5 * sin(π * (N-166)/80)
        elseif 247 ≤ N ≤ 366
            Eqt = 16.4 * sin(π * (N-247)/113)  
        end

        # Calculate the Solar Time (T_solar)
        T_solar = LT + Eqt/60 + (4*(standard_meridian_longitude - longitude))/60 # [DEGREES]
        solartimes[i] = T_solar
    end
    df[!, :SolarTime] = solartimes
    return df 
end    

function calculate_irradiance_new(solar_time, pyranometer, surface_azimuth, tilt_angle, k_t, local_time)

    # local_time = check_daylightsaving(local_time)
    N = dayofyear(local_time)
    T_solar = solar_time
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

    # Surface tilt β
    β = deg2rad(tilt_angle) # [RADIANS]
    # Surface azimuth γ
    γ = deg2rad(surface_azimuth) # [RADIANS]
    # Angle of Incidence (θ)
    θ = acos(sin(α) * cos(β) + cos(α) * sin(β) * cos(γ - Ψ)) # [RADIANS]

    # Direct Normal Irradiance 
    begin
        # Calculate diffuse fraction k_d by piecewise interpolation
        if k_t ≤ 0.22
            k_d = 1 - 0.09 * k_t
        elseif 0.22 < k_t ≤ 0.8
            k_d = 0.9511 - 0.1604 * k_t + 4.388 * k_t^2 - 16.638 * k_t^3 + 12.336 * k_t^4
        elseif k_t > 0.8
            k_d = 0.165
        end
    end

    # Calculate diffuse horizontal irradiance DHI
    DHI = pyranometer * k_d # [W/m^2]

    # Tilt angle of the pyranometer
    θz = deg2rad(0) # [RADIANS] 
    
    # Calculate direct normal irradiace DNI
    DNI = min((pyranometer - DHI)/cos(π/2 - α - θz), pyranometer) # [W/m^2]

    # If sun is not above horizon, all measured irradiance on the pyranometer must come from DHI
    if α < 0 
        DNI = 0
        DHI = pyranometer
    end

    # Beam Irradiance I_b
    # no direct beam irradiance if the angle of incidence is negative
    if cos(θ) < 0
        I_b = 0
    else
        I_b = DNI * cos(θ) # [W/m^2]
    end
    
    # GHI = DHI + DNI*(cos(π/2 - α))
    GHI = pyranometer

    # Diffuse irradiance I_d
    I_d = max(DHI * (1+cos(β))/2 + 0.5*GHI*(0.012*(π/2 - α)-0.04)*(1-cos(β)),0) # [W/m^2]
    
    # Reflected irradiance I_r
    albedo = 0.18 # albedo of grass
    I_r = max(GHI * albedo * (1 - cos(β)) / 2, 0) # [W/m^2]
    
    # Total irradiance I_total
    I_total = I_d + I_r + I_b # [W/m^2]
    return I_total
end

function calculate_solarheatgain(solar_time, GHI, k_t, local_time)
    
    Berg_tilt = 15 # [DEGREES] The Berg structure is tilted 15 degrees east of true south
    
    GHI = max(GHI, 0)

    I_t = calculate_irradiance_new(solar_time, GHI, 0, 0, k_t, local_time) # [W/m^2] top side

    I_e = calculate_irradiance_new(solar_time, GHI, 90-Berg_tilt, 90, k_t, local_time) # [W/m^2] east side

    I_s = calculate_irradiance_new(solar_time, GHI, 180-Berg_tilt, 90, k_t, local_time) # [W/m^2] south side

    I_w = calculate_irradiance_new(solar_time, GHI, 270-Berg_tilt, 90, k_t, local_time) # [W/m^2] west side

    I_n = calculate_irradiance_new(solar_time, GHI, 0-Berg_tilt, 90, k_t, local_time) # [W/m^2] north side

    # Calculate total solar heat gain from all 5 sides (3.412[BTU/hr/W], 10.764[ft^2/m^2])
    # TotalSolarHeatGain = 3.412*(I_t * L_wall * L_wall + (I_e + I_s + I_w + I_n) * L_wall * H_wall)/10.764 # [BTU/hr]
    # stick to the units [W]
    TotalSolarHeatGain = I_t * L_wall * L_wall + (I_e + I_s + I_w + I_n) * L_wall * H_wall # [W]

    return TotalSolarHeatGain
end

function estimate_cloudcover(Ta, Td, rh)
    # Calculate the temperature difference in Celcius
    temp_dew_diff = Ta - Td # [°C]
    
    # Using arbitrary thresholds to estimate cloud cover based on the temperature difference and relative humidity
    if temp_dew_diff < 2 && rh > 90  # 3.6°F approximates 2°C
        return 0.95  # [1] 90-100% cloud cover
    elseif temp_dew_diff < 4 && rh > 80  # 7.2°F approximates 4°C
        return 0.8 # [1] 70-89% cloud cover
    elseif temp_dew_diff < 6 && rh > 70  # 10.8°F approximates 6°C
        return 0.5  # [1] 30-69% cloud cover
    else
        return 0.15  # [1] 0-29% cloud cover
    end
end    

function estimate_Ts(Ti, Ta, solartime, day_status, GHI)
    if day_status == 0
        Ts = (Ta + Ti)/2
        return [Ts, Ts, Ts, Ts, Ts] # [°C]
    end

    # East side, peaks at 10 am
    Ts_E = max(Ti, Ti * ((GHI+1000)/2000) * (-13/7) * cos(2*pi*(solartime+2)/24)) # [°C]
    
    # South side, peaks at 2 pm
    Ts_S = max(Ti, Ti * ((GHI+1000)/2000) * (-13/7) * cos(2*pi*(solartime-2)/24)) # [°C]
    
    # West side, peaks at 4:30 pm
    Ts_W = max(Ti, Ti * ((GHI+1000)/2000) * (-13/7) * cos(2*pi*(solartime-4.5)/24)) # [°C]
    
    # North side, relatively cool
    Ts_N = max(Ti, Ti * (-8/7) * ((GHI+1000)/2000) * cos(2*pi*(solartime-48)/72)) # [°C]

    # Top side, not as hot, peaks at 2 pm
    Ts_R = max(Ti, Ti * ((GHI+3000)/4000) * (-12/7) * cos(2*pi*(solartime-52)/72)) # [°C]

    return [Ts_E, Ts_S, Ts_W, Ts_N, Ts_R] # [°C]
end

function Q_radiativeCooling(Ti, Ta, RH, GHI, solartime, day_status)

    # Calculate dew point temperature in Celcius
    temp_dewpoint = calculate_dew_point(Ta, RH) # [°C]

    # Emissivity of sky based on dew point temperature in Celcius
    e_sky = 0.741 + 0.0062 * temp_dewpoint # [1] 
    
    # Radiative Cooling using Stefan–Boltzmann Law
    sigma = 5.67*10^-8 # [W/(m^2*K^4)] Stefan–Boltzmann constant

    # Total exposed area subject to radiative cooling by side
    exposed_area = [L_wall*H_wall, L_wall*H_wall, L_wall*H_wall, L_wall*H_wall, L_wall*L_wall] # [m^2] 
    total_exposed_area = sum(exposed_area) # [m^2] 
    
    # Surface temperature by side 
    Ts_list = estimate_Ts(Ti, Ta, solartime, day_status, GHI) # [°C]
    Ts = Ts_list.+273.15 # [°K]

    T_sky = Ta * e_sky ^ 0.25 # [°C] # cite Genevieve paper

    # Calculate the net heat loss due to radiation which is the difference in radiation from Berg to ambient and radiation from ambient to Berg.
    # Q_radcool = 3.412 * sigma * (e_Berg * sum((Ts.^4).*exposed_area) - e_sky * total_exposed_area *(T_sky+273.15)^4) # [BTU/hr]
    # Stick to the unit [W]
    Q_radcool = sigma * (e_Berg * sum((Ts.^4).*exposed_area) - e_sky * total_exposed_area *(T_sky+273.15)^4) # [W]

    return Q_radcool
end    

function IsDay(solar_time, local_time)
    # Calculate the Solar Time (T_solar)
    T_solar = solar_time # [DEGREES]
    N = dayofyear(local_time)
    # Solar Declination (δ)
    δ = 23.45 * sin(deg2rad((360/365) * (284 + N))) # [DEGREES]

    # Solar Time Angle (ω)
    ω = 15 * (T_solar - 12) # [DEGREES]

    # Latitude (φ)
    φ = deg2rad(latitude) # [RADIANS]

    # Solar Elevation Angle (α)
    α = asin(sin(φ) * sin(deg2rad(δ)) + cos(φ) * cos(deg2rad(δ)) * cos(deg2rad(ω))) # [RADIANS]

    day_status = 0
    
    if α>0
        day_status = 1
    end
    
    return day_status
end    

# Function to normalize a single row with custom bounds
function normalize_row(row, cols, bounds)
    normalized_row = deepcopy(row)
    for col in cols
        if haskey(bounds, col)
            lower_bound, upper_bound = bounds[col]
            normalized_row[!, col] = (row[!, col] .- lower_bound) ./ (upper_bound - lower_bound)
        end
    end
    return normalized_row
end

# Function to find the nearest centroid
function find_nearest_centroid(normalized_row, centroids, cols)
    # Extract the relevant columns from the normalized_row and flatten the vector
    row_vector = vcat([normalized_row[!, col] for col in cols]...)
    
    # Calculate distances to each centroid
    distances = [norm(row_vector .- centroids[:, i]) for i in 1:size(centroids, 2)]
    
    # Return the index of the nearest centroid
    return argmin(distances)
end