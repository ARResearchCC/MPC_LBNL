# Model Predictive Control: LBNL Simulation. Version 3.0
# Complete Date: 06/18/2024
# ACE Microgrid Project --- CEE Atmosphere and Energy --- Stanford University

# This code is an intellectual property of Yuanbei "Fred" Fan.
# Dear user, this is Fred. If there are any questions with regards to the code, please do not feel hesistant to reach out to me at yf1098@stanford.edu.

############ Initialize Tools ############
begin
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

    include("PassiveModel.jl")
    include("MPCAlgorithms_LBNL.jl")
    ENV["PYTHON"] = "C:\\Users\\Fred\\anaconda3\\python.exe"
    ENV["GRB_LICENSE_FILE"] = "C:\\Users\\Fred\\.julia\\environments\\v1.7\\gurobi.lic" # Fred's license
    # ENV["GRB_LICENSE_FILE"] = "/Library/gurobi1100/gurobi.lic" # Andreas' license
    
    pushfirst!(PyVector(pyimport("sys")."path"), "C:\\Users\\Fred\\Desktop\\PyFMI")

    # Import the script as a module
    FMU = pyimport("FMU_Simulation_2")
end

############ Program Preparations ############
begin
    # Update automatically the date when this program is run.
    today_date = today()

    # Please update information of this program to automatically update the code name.
    version = 3.0

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
    # file_name = "C:\\Users\\Fred\\Desktop\\Microgrid_Project\\Data\\NREL_NSRDB_HalfMoonBay5_5Min\\103920_36.66_-121.01_2022.csv" # 5 mins data
    file_name = "C:\\Users\\Fred\\Desktop\\Microgrid_Project\\Data\\NREL_NSRDB_HalfMoonBay25_60Min_New\\137344_37.49_-122.42_2022.csv" # hourly data

    # Check if the file exists before trying to read it
    if isfile(file_name)
        weather = DataFrame(CSV.File(file_name, header=3))
    else
        println("File $file_name not found.")
    end
end

############ Declare Parameters ############
begin  
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
        # SetPointT_Low = 68 # [°F] Lower Bound of Set Point Temperature
        SetPointT_Low = 65 # [°F] Lower Bound of Set Point Temperature (To make sure the model doesn't run into infeasible situations, we might need to use this lower bound)
        SetPointT_High = 75 # [°F] Upper Bound of Set Point Temperature

        # Relative Humidity standard (currently not used due to lack of humidity data)
        # SetPointW_Low = 0.4 # LEED Regulation
        # SetPointW_High = 0.6 # LEED Regulation

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
        # TotalPersonHeat = 300 # [BTU/hr/PPL] Dan's suggestion
    end
    # Equipments Parameters
    begin
        # Device Capacity Parameters
        PVSize = 10; # [kW] PV DC Power Capacity
        BatterySize = 15; # [kWh] Battery Energy Capacity
        InverterSize = 15 # [kW] Max Continuous AC Output Power
        PCM_H_Size = 10 * 3412.14; # [BTU] PCM Heating Storage Energy Capacity
        PCM_C_Size = 1 * 3412.14; # [BTU] PCM Cooling Storage Energy Capacity

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
        BatteryLoss = 0.00001 # [/hr] Battery Leakage Rate
        MaxDischarge = 0.8 # [1] Max Depth of Discharge in Battery (80%)
        η = 0.98 # [1] Battery Inverter Efficiency

        # LBNL System Parameters
        # Heat Pump heating and cooling capacity
        Cap_HP_C = 3.85 # [kW] Max input electrical power for heat pump (cooling)
        Cap_HP_H = 4 # [kW] Max input electrical power for heat pump (heating)
        
        # Heat Pump heating and cooling coefficients of performance
        
        # HP_a = 6.08
        # HP_b = -0.0941
        # HP_c = 0.000464
        
        COP_H = 4 # COP of heating (currently constant)
        COP_C = 3.8 # COP of cooling (currently constant)

        C_HP_OP = 0.02 * Cap_HP_H # [$/(kW*YR)] Operational Cost of Heat Pump

        # Thermal Energy Storage - Phase Changing Materials
        Cap_PCM_H_Charging = Cap_HP_H * 4 * 3412.14 # [BTU/hr] Max Charging Power Capacity of PCM Heating Storage
        Cap_PCM_H_Discharging = Cap_PCM_H_Charging # [BTU/hr] Max Discharging Power Capacity of PCM Heating Storage
        Cap_PCM_C_Charging = Cap_HP_C * 4 * 3412.14 # [BTU/hr] Max Charging Power Capacity of PCM Cooling Storage
        Cap_PCM_C_Discharging = Cap_PCM_C_Charging # [BTU/hr] Max Discharging Power Capacity of PCM Cooling Storage

        η_PCM_H = 0.99 # [1] PCM Heating Storage Efficiency
        η_PCM_C = 0.99 # [1] PCM Cooling Storage Efficiency

        C_PCM_H_OP = 0.02 * PCM_H_Size # [$/(kWh*YR)] Operational Cost of PCM Heating Storage
        C_PCM_C_OP = 0.02 * PCM_C_Size # [$/(kWh*YR)] Operational Cost of PCM Cooling Storage

        # Standard Operating Power of Heat Pump and PCM Thermal Storages 
        P_H2HP = Cap_HP_H # [kW] default constant electrical power consumption for heat pump (heating)
        P_H2C = Cap_HP_C # [kW] default constant electrical power consumption for heat pump (cooling)
        P_HP2PCM_H = Cap_PCM_H_Charging # [BTU/hr] default constant heat charging rate of PCM Heating Storage
        P_C2PCM_C = Cap_PCM_C_Charging # [BTU/hr] default constant heat charging rate of PCM Cooling Storage
        P_PCM_H2H = Cap_PCM_H_Discharging # [BTU/hr] default constant heat discharging rate of PCM Heating Storage
        P_PCM_C2H = Cap_PCM_C_Discharging # [BTU/hr] default constant heat discharging rate of PCM Cooling Storage
        
        # Big M Method
        M = 10000
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
    # Make SimpleSchedule
    begin
        SimpleSchedule = zeros(8760, 3);
        for i = 0:364
            for j = 10:17
                SimpleSchedule[i*24 + j, 1] = 1
                SimpleSchedule[i*24 + j, 2] = 1
                SimpleSchedule[i*24 + j, 3] = 1
            end
        end
    end
    # Make ComplexSchedule
    begin
        ComplexSchedule = zeros(8760, 3);
        for i = 0:364
            for j = 1:6
                ComplexSchedule[i*24 + j, 1] = 0
                ComplexSchedule[i*24 + j, 2] = 0
                ComplexSchedule[i*24 + j, 3] = 1
            end    
            for j = 7:8
                ComplexSchedule[i*24 + j, 1] = 1
                ComplexSchedule[i*24 + j, 2] = 0
                ComplexSchedule[i*24 + j, 3] = 1
            end  
            for j = 9:19
                ComplexSchedule[i*24 + j, 1] = 1
                ComplexSchedule[i*24 + j, 2] = 0
                ComplexSchedule[i*24 + j, 3] = 0.5
            end  
            for j = 20:20
                ComplexSchedule[i*24 + j, 1] = 1
                ComplexSchedule[i*24 + j, 2] = 1
                ComplexSchedule[i*24 + j, 3] = 0.5
            end  
            for j = 21:22
                ComplexSchedule[i*24 + j, 1] = 1
                ComplexSchedule[i*24 + j, 2] = 1
                ComplexSchedule[i*24 + j, 3] = 1
            end  
            for j = 23:24
                ComplexSchedule[i*24 + j, 1] = 0
                ComplexSchedule[i*24 + j, 2] = 1
                ComplexSchedule[i*24 + j, 3] = 1
            end  
        end
    end
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
        for i = 1:8760
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

# Time Horizon Parameters
begin
    TimeStart = 1;
    TimeEnd = 8760;
    starttime = weather.datetime[TimeStart];
    f_run = 1 # run frequency [hour/time]
    Opt_Horizon = 168 # [hr]
    stepsize = 60*(1/f_run) # [min]
    δt = stepsize/60 # [hr] Declare stepzize for the optimization program
    # NumRun = (TimeEnd-TimeStart+1) - Opt_Horizon + 1; # This is the max number of runs allowed for a full year of data
    NumRun = 5
end   

############ Initialize Space ############
begin
    # Initialize M_States
    global current_m_states = DataFrame(Indoor_Temp = [294.15], PCM_Hot_Temp = [273.15+48.9], PCM_Cold_Temp = [273.15+10], 
    PCM_Hot_SOC = [1], PCM_Cold_SOC = [0.983103]);
    
    # Initialize J_States: [PV_generation, E_0, InStorageBattery]
    global current_j_states = [0, 0, 0.5*BatterySize]
    
    Costs = zeros(NumRun);
    AccumulatedCosts = zeros(NumRun);
    # AllStates = zeros(NumRun, 4);

    global TotalCost = 0;

    # ActionMatrix = zeros(NumRun, 13);
end

############ Program Execution ############
runtime = @elapsed begin
    for i = TimeStart:f_run:NumRun
        
        # Define timepoints needed for FMU simulation
        timepoints = [Dates.value(weather.datetime[i]-starttime)/1000, Dates.value(weather.datetime[i+1]-starttime)/1000] # [s]
        
        # Call the MPC function
        currentcost, new_j_states, input_schedule = Optimize1(i, weather[i : i+Opt_Horizon-1, :], SimpleSchedule[i : i+Opt_Horizon-1, :], current_m_states, current_j_states);
        
        
        println()
        println("Iteration")
        println(i)
        println()
        println("Input schedule")
        println(input_schedule)
        println()
        println("timepoints")
        println(timepoints)
        println("Calling FMU")

        # Collect the M_States at the final timestep of last FMU iteration and feed them as the initial M_States of the current FMU iteration.
        msize = size(current_m_states)[1]
        Temp_PCM_C_1 = current_m_states[msize, "PCM_Cold_Temp"]  # [°K]
        Temp_PCM_H_1 = current_m_states[msize, "PCM_Hot_Temp"]  # [°K]
        Temp_Indoor_1 = current_m_states[msize, "Indoor_Temp"]  # [°K]
        m_states_updates = [Temp_PCM_H_1, Temp_PCM_C_1, Temp_Indoor_1]

        # Call the FMU function with the input (current iteration, operational commands, timepoints)
        res = FMU.fmu(i, m_states_updates, input_schedule, timepoints)

        # This is the dataframe output from the current FMU iteration with power consumption and soc updates.
        new_m_states = pandas_to_julia(res)
        
        # Updates the Julia and Modelica states.
        global current_m_states = new_m_states
        global current_j_states = new_j_states
        
        # Add the cost to total cost
        Costs[i] = currentcost;
        AccumulatedCosts[i] = sum(Costs[n] for n = 1:i)
        global TotalCost += currentcost;
    end
end

# Plot
trace1 = scatter(
    x = 1:NumRun,  
    y=AccumulatedCosts,
    name="Optimized Storage System"
)

p_ll = plot([trace1], Layout(title="Accumulated Loss of Load Over Time", xaxis_title="Hours", yaxis_title="Loss of Load (kWh)"))
display(p_ll)
save_plot(p_ll, folder_path, "Accumulated Loss of Load MPC $version)_$today_date", "png")

println("The runtime is $runtime seconds")

println("The total loss of load is $TotalCost kWh")
