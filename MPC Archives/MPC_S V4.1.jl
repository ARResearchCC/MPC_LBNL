# Model Predictive Control: LBNL Simulation. Version 4.1
# Complete Date: 08/04/2024
# ACE Microgrid Project --- CEE Atmosphere and Energy --- Stanford University

# This code is an intellectual property of Yuanbei "Fred" Fan.
# Dear user, this is Fred. If there are any questions with regards to the code, please do not feel hesistant to reach out to me at yf1098@stanford.edu.

############ Initialize Tools ############
begin
    using Pkg

    # Please unmute the next two lines if you are running this code for the first time and don't have all the packages installed.
    # Pkg.add("JuMP", "CSV", "DataFrames", "PlotlyJS", "Dates", "XLSX", "FileIO", "Base", "Random", "Statistics", "Gurobi", "PyCall")
    # Pkg.update("Gurobi")
    
    using JuMP
    using CSV
    using DataFrames
    using PlotlyJS
    using Dates
    using XLSX
    using FileIO
    using Base
    using Random
    using Statistics
    using Gurobi
    using PyCall

    # Code setup
    include("PassiveModel.jl")
    include("MPCAlgorithms_LBNL_New.jl")
    ENV["PYTHON"] = "C:\\Users\\Fred\\anaconda3\\python.exe"
    ENV["GRB_LICENSE_FILE"] = "C:\\Users\\Fred\\gurobi.lic" # Fred's license

    pushfirst!(PyVector(pyimport("sys")."path"), "C:\\Users\\Fred\\Desktop\\PyFMI")
    ############ Load FMU ############
    begin
        pyf = pyimport("pyfmi")
        FMU = pyimport("FMU_Simulation_4")

        # Define the path to the FMU file
        fmuName = "C:/Users/Fred/Desktop/PyFMI/Stanford_Hybrid_System.fmu"
        # Load the FMU model
        model = pyf.load_fmu(fmuName)
    end
end

############ Program Preparations ############
begin
    # Update automatically the date when this program is run.
    today_date = today()

    # Please update information of this program to automatically update the code name.
    version = 4.1

    code_name = "MPC_S_V$version._$today_date"

    # Create folder to later save data and plots
    begin
        # Define the path for the new folder
        folder_path = "C:\\Users\\Fred\\Desktop\\Microgrid_Project\\Test_Results\\$code_name"

        # Use mkpath() to create the folder if it doesn't exist
        mkpath(folder_path)
    end

    # Function to save a plot as a PNG file in the specified folder
    function save_plot(plot, path, filename, format="png")
        # Create the full file path with the specified filename and format
        full_path = joinpath(path, string(filename, ".", format))
        
        # Save the plot as an image in the desired format
        savefig(plot, full_path)
    end
end

############ Load Data ############
begin
    # Load the weather file. This weather file is a historical weather file with 5 minutes step size. This weather file is treated as the ground truth, but simulated forecast errors are added later.
    file_name = "C:\\Users\\Fred\\Desktop\\Microgrid_Project\\Data Cleaning\\updated_epw_weather_nonleap_5mins.csv"
    # file_name = "C:\\Users\\Fred\\Desktop\\Microgrid_Project\\Data\\HMB_data_5mins\\57357_37.46_-122.43_2022.csv"

    # Check if the file exists before trying to read it
    if isfile(file_name)
        # weather = DataFrame(CSV.File(file_name, header=3))
        weather = DataFrame(CSV.File(file_name))
    else
        println("File $file_name not found.")
    end
end

############ Declare Parameters ############
begin 
    # Time Horizon Parameters
    begin
        TimeStart = 1;
        f_run = 12 # run frequency [times/hour]
        stepsize = 60*(1/f_run) # [min]
        TimeEnd = 8760 * f_run;
        # Opt_Horizon = 6 * f_run # [intervals]
        # NumRun = (TimeEnd-TimeStart+1) - Opt_Horizon + 1; # This is the max number of runs allowed for a full year of data
        NumRun = 3

        # Define the sampling intervals
        stepsizes = [5, 30, 60, 120] # minutes
        steps = Int[]
        append!(steps, fill(stepsizes[1], 6))
        # append!(steps, fill(stepsizes[2], 5))
        # append!(steps, fill(stepsizes[3], 21))
        # append!(steps, fill(stepsizes[4], 48))
    end    
    # Berg Envelope Parameters
    begin
        L_wall = 20 # [ft] Length of the Wall
        H_wall = 8 # [ft] Height of the Wall
        R_wall = 15 # [ft^2·°F·hr/BTU] Thermal Resistance of the Wall
        R_floor = 17 # [ft^2·°F·hr/BTU] Thermal Resistance of the Floor
        H_Ceiling = 80/12 # [ft] Interior Height of the Ceiling (80 inch) 
        Area = L_wall * L_wall # [ft^2] Floor/Ceiling Area 
        Berg_tilt = 15 # [DEGREES] The Berg structure is tilted 15 degrees east of true south
        UA = L_wall*H_wall*4*(1/R_wall) + Area*2*(1/R_floor) # [BTU/(hr⋅°F)]
        TC = 581.7 # [BTU/°F] Total Thermal Capacitance of the Structure - calibrated using passive model 4.6
        Volume = 2330 # [ft^3]
        #=
        UA = 89, which is consistent with given documents. Note that this is excluding infiltration, since infiltration is considered 
        as a heat transfer q_Infiltration instead of extra UA.
        =#
    end
    # Radiation Parameters
    begin
        e_Berg = 0.75 # currently not used
        SHGC = 0.0215 # calibrated using passive model 4.6
        RCC = 0.0408 # calibrated using passive model 4.6
    end
    # Environmental Parameters
    begin
        wind_height = 9.144 # [m] The height above ground at which wind speed is measured. The PVWatts default is 9.144 m.
        albedo = 0.18 # [1] albedo of grass
        n_air = 1 # [1] Index of reflaction of air
        n_AR = 1.3 # [1] Index of reflaction of AR coating
        n_glass = 1.526 # [1] Index of reflaction of glass as in PVWatt model for standard modules
        standard_meridian_longitude = -120; # [DEGREES] Longitude of standard meridian
        longitude = -122.4286 # [DEGREES] Longitude of site
        latitude = 37.4636 # [DEGREES] Latitude of site
    end
    # Infiltration Parameters
    begin
        # Stack coefficient for building(house height = 1 story)
        Cs = 0.000145 # [(L/s)^2/(cm^4*K)]
        # Wind coefficient for building(house height =1 story, shelter class=1(no obstructions or local shielding))
        Cw = 0.000319 # [(L/s)^2/(cm^4*(m/s)^2)] (between 0.000319 and 0.000246 (house height =1 story, shelter class=2(typical shelter for an insolated rural house)))
        
        # Effective leakage area measured during Blower Door Test 
        ELA =  38.3 # [in^2] (between 38.3(Dan's interpolation) and 47.1(Jessie's interpolation))
        Al = ELA * 6.4516 # [cm^2]

        T_d = 273.15 - 1.67 # [°K] 99% Design Temperature for Half Moon Bay
        T_indoor_constant = 22 # [°C] Constant Indoor Temperature (for the simplicity of a linear model)
    end
    # Daylight Saving
    begin
        dls_start = DateTime(2024, 3, 10, 2)
        dls_end = DateTime(2024, 11, 3, 0)
        timezone = "America/Los_Angeles"  # Specify local timezone
    end
    # Interior Standard
    begin
        # Temperature standard
        # Temperature setpoints control the direct heating or cooling of the space through heat pump or PCM discharge.
        zone_temp_heating_setpoint = 20 * (9/5) + 32 # [°F] (20 °C)
        zone_temp_cooling_setpoint = 22 * (9/5) + 32 # [°F] (22 °C)
        zone_temp_setpoint_delta = 0.5 * (9/5) # [°F] (0.5 °C)

        # Ventilation
        Ventilation = 15 # [CFM/PPL] Building Standard
        Ra = 0.06 # [CFM/ft^2] Area Outdoor Air Rate; Source:CEE226E Week 5 - Energy Modeling Questions - Slide 9
        Rp = 5 # [CFM/PPL] People Outdoor Air Rate; Source:CEE226E Week 5 - Energy Modeling Questions - Slide 9   
    end
    # Lighting, Plugs, and Occupancy
    begin
        PeakLighting = 0.1 # [KW] Dan's suggestion
        PeakPlugLoad = 0.1 # [KW] Computer = 40W, Phone = 10W, 2*(40 + 10) = 100 [W] = 0.1 [kW]
        MaxOccupancy = 4 # [PPL] 4-MAN Office Room Plan
        PersonLatentHeat = 200; # [BTU/hr/PPL]  CEE226E Slide
        PersonSensibleHeat = 300; # [BTU/hr/PPL]  CEE226E Slide
        TotalPersonHeat = PersonSensibleHeat + PersonLatentHeat # [BTU/hr/PPL]
    end
    # Equipments Parameters
    begin
        # Device Capacity Parameters
        PVSize = 20; # [kW] PV DC Power Capacity
        BatterySize = 15; # [kWh] Battery Energy Capacity
        InverterSize = 15 # [kW] Max Continuous AC Output Power
        PCM_H_Size = 11 * 3412.14; # [BTU] PCM Heating Storage Energy Capacity
        PCM_C_Size = 7 * 3412.14; # [BTU] PCM Cooling Storage Energy Capacity

        # Solar PV Parameters
        noct_installed = 45 # [°C] The “installed” nominal operating cell temperature. PVWatts assumes this value to be 45 C for rack-mounted arrays and 49 C for roof mount systems with restricted air flow around the module.
        module_height = 5 # [m] The height above ground of the center of the module. The PVWatts default is 5.0.
        module_width = 0.31579 # [m] Module width. The default value of 0.31579 meters in combination with the default module_length gives a hydraulic diameter of 0.5.
        module_length = 1.2 # [m] Module length. The default value of 1.2 meters in combination with the default module_width gives a hydraulic diameter of 0.5.
        module_emissivity = 0.84 # [1] The effectiveness of the module at radiating thermal energy.
        module_absorption = 0.83 # [1] The fraction of incident irradiance that is converted to thermal energy in the module. 
        module_surface_tilt = 27 # [DEGREES] Module tilt from horizontal. If not provided, the default value of 30 degrees is used.
        Γ_t = -0.47/100 # [1/°C] Temperature coefficient for standard module
        T_ref = 25 # [°C] Reference cell temperature
        η_PV = 0.86 # [1] PV DC efficiency after system loss
        P_dc0 = 1 # [kW/kW] Rated capacity at standard conditions
        η_PVIV = 0.94 # [1] PV(DC) to Home(AC) inverter efficiency
        
        # Battery parameters
        BatteryLoss = 0.00001 # [/hr] Battery Leakage Rate ??
        MaxDischarge = 0.8 # [1] Max Depth of Discharge in Battery (80%)
        η = 0.98 # [1] Battery Inverter Efficiency

        # LBNL System Parameters
        
        # Heat Pump heating and cooling coefficients of performance
        
        # HP_a = 6.08
        # HP_b = -0.0941
        # HP_c = 0.000464
        
        COP_H = 4 # COP of heating (currently constant)
        COP_C = 3.8 # COP of cooling (currently constant)

        # C_HP_OP = 0.02 * Cap_HP_H # [$/(kW*YR)] Operational Cost of Heat Pump

        # η_PCM_H = 1/24 # [/hr] PCM Heating Storage Leakage Rate
        # η_PCM_C = 1/24 # [/hr] PCM Cooling Storage Leakage Rate

        C_PCM_H_OP = 0.02 * PCM_H_Size # [$/(kWh*YR)] Operational Cost of PCM Heating Storage
        C_PCM_C_OP = 0.02 * PCM_C_Size # [$/(kWh*YR)] Operational Cost of PCM Cooling Storage

        # Standard Operating Power of Heat Pump and PCM Thermal Storages 
        HP_power_H = 2 # [kW] default constant electrical power consumption for heat pump (heating)
        HP_power_C = 2 # [kW] default constant electrical power consumption for heat pump (cooling)
        PCM_H_discharge_rate = 2 * 3412.14 # [BTU/hr] default constant heat discharging rate of PCM Heating Storage
        PCM_C_discharge_rate = 2 * 3412.14 # [BTU/hr] default constant heat discharging rate of PCM Cooling Storage
        
        # Standard Operating Power of Control System, Air Handling Unit (AHU), and Pump
        P_Controls = 0.02 # [kW]
        P_AHU = 0 # [kW]
        P_Pumps = 0.04 # [kW]

        # Big M Method
        M = 200
    end
end

############ Generate Schedules ############
# We should add bigger and maybe more variable electrical load profiles
begin
    # Instructions
    begin
        # The are two schedules, format: Schedules[Lighting, Plugs, Occupancy]
        # All variables are unitless(0.5 = 50%). They are then used to multiply peak lighting/occupancy/plugs heat gain in the code.

        # SimpleSchedule:
        # About lighting schedule: full lighting intensity from 6 am to 10 pm daily
        # About plugs schedule: full plugs intensity from 6 am to 10 pm daily
        # About occupancy schedule: full occupancy from 6 am to 10 pm daily
    
        # ComplexSchedule:
        # About lighting schedule: full lighting intensity from 6 am to 10 pm daily
        # About plugs schedule: equipments (PC and phone) are charged from 7 pm to 12 am daily (5 hours)
        # About occupancy schedule: 50% occupany from 8 am to 8 pm, 100% occupancy during night time.
    end
    function generate_schedules(f_run::Int)
        # Calculate the total number of intervals in a year
        total_intervals = TimeEnd
        
        # Initialize schedules
        SimpleSchedule = zeros(total_intervals, 3)
        ComplexSchedule = zeros(total_intervals, 3)
        
        # Fill SimpleSchedule
        for i in 0:364
            for j in 10*f_run:(17*f_run-1)
                SimpleSchedule[i*24*f_run + j + 1, :] .= 1
            end
        end
        
        # Fill ComplexSchedule
        for i in 0:364
            for j in 1*f_run:(6*f_run-1)
                ComplexSchedule[i*24*f_run + j + 1, 1:2] .= 0
                ComplexSchedule[i*24*f_run + j + 1, 3] = 1
            end
            for j in 7*f_run:(8*f_run-1)
                ComplexSchedule[i*24*f_run + j + 1, 1] = 1
                ComplexSchedule[i*24*f_run + j + 1, 2] = 0
                ComplexSchedule[i*24*f_run + j + 1, 3] = 1
            end
            for j in 9*f_run:(19*f_run-1)
                ComplexSchedule[i*24*f_run + j + 1, 1] = 1
                ComplexSchedule[i*24*f_run + j + 1, 2] = 0
                ComplexSchedule[i*24*f_run + j + 1, 3] = 0.5
            end
            for j in 20*f_run:(20*f_run+f_run-1)
                ComplexSchedule[i*24*f_run + j + 1, 1:2] .= 1
                ComplexSchedule[i*24*f_run + j + 1, 3] = 0.5
            end
            for j in 21*f_run:(22*f_run-1)
                ComplexSchedule[i*24*f_run + j + 1, :] .= 1
            end
            for j in 23*f_run:(24*f_run-1)
                ComplexSchedule[i*24*f_run + j + 1, 1] = 0
                ComplexSchedule[i*24*f_run + j + 1, 2] = 1
                ComplexSchedule[i*24*f_run + j + 1, 3] = 1
            end
        end
    
        return SimpleSchedule, ComplexSchedule
    end
    
    SimpleSchedule, ComplexSchedule = generate_schedules(f_run)
    
end

############ Data Preparations ############
begin
    # Corrected function to create DateTime from row
    begin
        # Apply function to each row to create a new column with DateTime objects
        weather.datetime = [create_datetime_from_row(row) for row in eachrow(weather)] 
    end
    
    # Use PyCall to format data into Pandas
    begin
        # Import necessary Python libraries
        pd = pyimport("pandas")
    
        # Convert Julia DataFrame to Pandas DataFrame
        columns = names(weather)
        data = [getproperty(weather, col) for col in columns]
        pandas_df = pd.DataFrame(Dict(zip(columns, data)))
    
        # Ensure the datetime column is in the correct format and set as the index
        pandas_df["datetime"] = pd.to_datetime(pandas_df["datetime"])
        pandas_df = pandas_df.set_index("datetime")
    end
    
    # Calculate PV module cell temperature with PVLib [°C]
    begin
        # Proceed with pvlib operations
        pvlib = pyimport("pvlib")
    
        # Prepare data input for pvlib cell temperature calculation
        ambient_temperature = pandas_df["Temperature"] # [°C]
        wind_speed = pandas_df["Wind Speed"] # [m/s]
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
    
    # Calculate PV output [kW DC Output/kW Capacity]
    begin
        # Calculate PV output [kW DC Output/kW Capacity]  
        PV_output = [];
        for i = 1:TimeEnd
            pv_o = calculate_PV(weather[i, :datetime], weather[i, :DNI], weather[i, :DHI], weather[i, :GHI], 180-Berg_tilt, module_surface_tilt, weather[i, :celltemp])
            push!(PV_output, pv_o)
        end
        
        # Store PV output into weather file
        weather.PV = PV_output
    end

    # Add solar time column to the weather file
    weather = add_solar_time(weather)
    
    # Convert [°C] to [°F]
    weather.Temperature .= (weather.Temperature .* (9/5)) .+ 32 # [°F]
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

function p2e(M_States)
    
    # Length of the M_States dataframe
    msize = size(M_States)[1];
    
    time_intervals = zeros(msize-1); # [hr]

    for i = 1:msize-1
        time_intervals[i] = (M_States[i+1, :"Time"] - M_States[i, :"Time"])/3600 # [hr]
    end

    print(time_intervals)

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

############ Initialize Space ############
begin
    # Initialize M_States to start the FMU simulation and Julia optimization
    #=
    These are the current initial conditions:
        PCM Hot Initial Temperature = 48 [°C]
        PCM Cold Initial Temperature = 11 [°C]
        Initial Indoor Temperature = 22 [°C]
        PCM Hot Capacity = PCM_H_Size/3412.14 = 11 [kWh]
        PCM Cold Capacity = PCM_C_Size/3412.14 = 7 [kWh]
        Zone Temperature Cooling Setpoint = 22 [°C]
        Zone Temperature Heating Setpoint = 20 [°C]
        Zone Temperature Delta Setpoint = 0.5 [°C]
        PCM Cold Initial SOC = 0.5
        PCM Hot Initial SOC = 0.5
    =#

    global current_m_states = DataFrame(PCM_Hot_Temp = [273.15+48], PCM_Cold_Temp = [273.15+11], Indoor_Temp = [273.15+22], 
    PCM_Hot_Cap = [PCM_H_Size/3412.14], PCM_Cold_Cap = [PCM_C_Size/3412.14],  
    T_Cooling = [273.15+(zone_temp_cooling_setpoint-32)*(5/9)], T_Heating = [273.15+(zone_temp_heating_setpoint-32)*(5/9)], 
    T_Delta = [(zone_temp_setpoint_delta)*(5/9)], PCM_Cold_SOC = [0.5], PCM_Hot_SOC = [0.5]);
    
    # This is the same as current_m_states at i=1, but it works with the python input for the FMU.
    initialstates = [];
    for i = 1:size(current_m_states)[2]
        push!(initialstates, current_m_states[1, i])
    end

    # Initialize J_States: [PV_generation, E_0, InStorageBattery]
    #=
    These are the current initial conditions:
        PV Generation at last timestep = 0 [kW]
        Actual Power consumption at last timestep = 0 [kW]
        Battery Initial SOC = 0.5
    =#
    global current_j_states = [0, 0, 0.5*BatterySize]
    
    Costs = zeros(NumRun);
    AccumulatedCosts = zeros(NumRun);
    fmu_temp = zeros(NumRun);
    fmu_indoor_temp = zeros(NumRun);
    Solve_failure = zeros(NumRun);
    AllStates = zeros(NumRun, 4);

    Heat_pump_modes_j = zeros(NumRun);
    HP_H_Usage_j = zeros(NumRun)
    
    fmu_HP_Mode = zeros(NumRun)
    fmu_HP_power = zeros(NumRun)
    fmu_Useful_Thermal= zeros(NumRun)
    fmu_HT_PCM_C = zeros(NumRun)
    fmu_HT_PCM_H = zeros(NumRun)
    fmu_system_Mode = zeros(NumRun)
    fmu_FCU = zeros(NumRun)
    fmu_fan = zeros(NumRun)
    fmu_Pump1 = zeros(NumRun)
    fmu_Pump2 = zeros(NumRun)

    global TotalCost = 0;

    Actions = zeros(NumRun);
end

############ Program Execution ############
runtime = @elapsed begin
    starttime = weather.datetime[TimeStart];
    
    initial_schedule = [0, 0]; # Do nothing while the model is solving the first iteration
    initial_timepoints = [0, 60*stepsize]
    res = FMU.fmu(1, initialstates, initial_schedule, initial_timepoints, model)
    initial_m_states = pandas_to_julia(res)
    
    Actions[1] = 0;   
    AllStates[1, 1] = current_j_states[3]
    AllStates[1, 2] = initial_m_states[end, "PCM_Hot_SOC"]
    AllStates[1, 3] = initial_m_states[end, "PCM_Cold_SOC"]
    
    fmu_temp[1] = initial_m_states[1, "Ambient Temp"]
    fmu_indoor_temp[1] = initial_m_states[end, "Indoor_Temp"]
    fmu_HP_Mode[1] = initial_m_states[1, "HP Mode"]
    fmu_system_Mode[1] = initial_m_states[1, "Exact System Mode"]
    
    fmu_fan[1] = p2e(initial_m_states)[1]
    fmu_FCU[1] = p2e(initial_m_states)[2]
    fmu_HP_power[1]= p2e(initial_m_states)[3]
    fmu_HT_PCM_C[1] = p2e(initial_m_states)[4]
    fmu_HT_PCM_H[1] = p2e(initial_m_states)[5]
    fmu_Pump1[1] = p2e(initial_m_states)[6]
    fmu_Pump2[1] = p2e(initial_m_states)[7]
    fmu_Useful_Thermal[1]= p2e(initial_m_states)[8]

    # Updates the Julia and Modelica states.
    global current_m_states = initial_m_states
    

    for i = TimeStart+1:1:NumRun
        
        println()
        println("Iteration: ", i)
        println()
        println("The current M states are:")
        println(current_m_states)

        # Define timepoints needed for FMU simulation
        timepoints = [Dates.value(weather.datetime[i]-starttime)/1000, Dates.value(weather.datetime[i+1]-starttime)/1000] # [s]
        # println(timepoints)

        input_data = passive_model(i, weather::DataFrame, ComplexSchedule, steps, stepsize)
        println("The input data from the passive model is:")
        println(input_data)
        # Call the MPC function
        currentcost, new_j_states, input_schedule, didnt_solve = OptimizeV4_1(i-1, steps, stepsize, input_data, current_m_states, current_j_states);

        println("Optimization Successful!")
        println()
        println("Input schedule: ", input_schedule)
        println()
        println("timepoints: ", timepoints)
        println()
        println("Calling FMU")

        # Call the FMU function with the input (current iteration, operational commands, timepoints, FMU model)
        res = FMU.fmu(i, initialstates, input_schedule, timepoints, model)
        println()
        println("The Modelica Feedbacks are:")
        println(res)
        # This is the dataframe output from the current FMU iteration with power consumption and SOC updates etc. (Modelica states)
        new_m_states = pandas_to_julia(res)
        
        # Updates the Julia and Modelica states.
        global current_m_states = new_m_states
        global current_j_states = new_j_states
        
        # Add the lost of load (cost) from the current iteration to total lost of load (cost)
        Costs[i] = currentcost;
        AccumulatedCosts[i] = sum(Costs[n] for n = 1:i)
        
        # Store return ambient temperature from fmu to verify epw file
        fmu_temp[i] = new_m_states[1, "Ambient Temp"]
        fmu_indoor_temp[i] = new_m_states[end, "Indoor_Temp"]
        fmu_HP_Mode[i] = new_m_states[1, "HP Mode"]
        fmu_system_Mode[i] = new_m_states[1, "Exact System Mode"]
        
        fmu_fan[i] = p2e(new_m_states)[1]
        fmu_FCU[i] = p2e(new_m_states)[2]
        fmu_HP_power[i]= p2e(new_m_states)[3]
        fmu_HT_PCM_C[i] = p2e(new_m_states)[4]
        fmu_HT_PCM_H[i] = p2e(new_m_states)[5]
        fmu_Pump1[i] = p2e(new_m_states)[6]
        fmu_Pump2[i] = p2e(new_m_states)[7]
        fmu_Useful_Thermal[i]= p2e(new_m_states)[8]

        # Is the current iteration of optimization infeasible?
        Solve_failure[i] = didnt_solve
        
        # Store the battery SOC, PCM_H SOC, PCM_C SOC, and 
        
        AllStates[i, 1] = new_j_states[3]
        AllStates[i, 2] = new_m_states[1, "PCM_Hot_SOC"]
        AllStates[i, 3] = new_m_states[1, "PCM_Cold_SOC"]

        # Store actions
        Actions[i] = input_schedule[1]
        global TotalCost += currentcost;
    end
end

#=
# Plot Accumulated loss of load
trace1 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=AccumulatedCosts,
    name="Optimized Storage System"
)

p_ll = plot([trace1], Layout(title="Accumulated Loss of Load Over Time", xaxis_title="Hours", yaxis_title="Loss of Load (kWh)"))
display(p_ll)
# save_plot(p_ll, folder_path, "Accumulated Loss of Load MPC $version)_$today_date", "png")
=#

#=
# Plot Temperature Anamoly
trace2 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=fmu_temp.-fahrenheit_to_kelvin(weather[1:NumRun, "Temperature"]),
    name="Temp Diff"
)

p_t = plot([trace2], Layout(title="Temp Difference", xaxis_title="Time Steps", yaxis_title="Temp Anamoly (K)"))
display(p_t)
=#

trace7 = scatter(
    x = 1:NumRun,  
    y=Solve_failure,
    name="Solve Failure"
)
p_f = plot([trace7], Layout(title="Solve Failure vs. Time", xaxis_title="Time Steps", yaxis_title=""))
display(p_f)

# Plot Indoor and Ambient Temperatures
trace8 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=fmu_indoor_temp.-273.15,
    name="Indoor Temp"
)

trace9 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=fmu_temp.-273.15,
    name="Ambient Temp"
)

p_ts = plot([trace8, trace9], Layout(title="Temperatures", xaxis_title="Time Steps", yaxis_title="Temperatures (C)"))
display(p_ts)

# Plot Julia vs. FMU actions
trace3 = scatter(
    x = 1:NumRun,  
    y=Actions,
    name="Julia Command"
)

trace12 = scatter(
    x = 1:NumRun,  
    y=fmu_system_Mode,
    name="FMU System Mode"
)

trace19 = scatter(
    x = 1:NumRun,  
    y=fmu_HP_Mode,
    name="FMU HP Mode"
)

p_a = plot([trace3, trace12, trace19], Layout(title="HP Mode, Julia Command vs. FMU Actions", xaxis_title="Time Steps", yaxis_title=""))
display(p_a)

# Plot SOCs
trace4 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=AllStates[:, 1]/BatterySize,
    name="Battery SOC"
)
trace5 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=AllStates[:, 2],
    name="PCM H SOC"
)
trace6 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=AllStates[:, 3],
    name="PCM C SOC"
)

p_soc = plot([trace4, trace5, trace6], Layout(title="SOC vs. Time", xaxis_title="Time Steps", yaxis_title="SOC percentage"))
display(p_soc)

trace10 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=fmu_fan,
    name="FMU Fan Energy Consumption"
)

trace11 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=fmu_FCU,
    name="FMU Fan Coil Unit Heat Delivered"
)

trace13 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=fmu_HP_power,
    name="fmu HP Energy Consumption"
)

trace14 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=fmu_HT_PCM_C,
    name="FMU Heat Transfer PCM C"
)

trace15 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=fmu_HT_PCM_H,
    name="FMU Heat Transfer PCM H"
)

trace16 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=fmu_Pump1,
    name="FMU Pump 1"
)

trace17 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=fmu_Pump2,
    name="FMU Pump 2"
)

trace18 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=fmu_Useful_Thermal,
    name="FMU Useful Thermal Delivered"
)

trace19 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=AllStates[:, 2].*(PCM_H_Size/3412.14),
    name="PCM H Remaining Energy"
)
trace20 = scatter(
    x = weather[1:NumRun, "datetime"],  
    y=AllStates[:, 3].*(PCM_C_Size/3412.14),
    name="PCM C Remaining Energy"
)

p_energy = plot([trace10, trace16, trace17], Layout(title="Device Energy Consumption every 5 minutes", xaxis_title="Time Steps", yaxis_title="Energy Consumption (kWh)"))
display(p_energy)

p_thermal = plot([trace18, trace11, trace13], Layout(title="Device Thermal Energy every 5 minutes", xaxis_title="Time Steps", yaxis_title="Energy (kWh)"))
display(p_thermal)

p_pcm = plot([trace14, trace15, trace19, trace20], Layout(title="PCM Plots", xaxis_title="Time Steps", yaxis_title="Energy (kWh)"))
display(p_pcm)


# Runtime per iteration (must be less than the MPC run frequency)
rt = runtime/(NumRun-TimeStart)
println("The runtime is $rt seconds")

println("The total loss of load is $TotalCost kWh")

Total_Inf = sum(Solve_failure)
Inf_pct = 100*Total_Inf/(NumRun-1)

println("There are a total of $Total_Inf infeasible iterations out of $(NumRun-1) iterations, which is $Inf_pct % of the total iterations.")