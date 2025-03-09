# Rule Based Model Predictive Control: LBNL Simulation. Version 1.0
# Complete Date: 03/05/2025
# ACE Microgrid Project --- CEE Atmosphere and Energy --- Stanford University

# This code is an intellectual property of Yuanbei "Fred" Fan.
# Dear user, this is Fred. If there are any questions with regards to the code, please do not feel hesistant to reach out to me at yf1098@stanford.edu.

############ Initialize Tools ############
begin
    using Pkg

    # Please unmute the next two lines if you are running this code for the first time and don't have all the packages installed.
    # Pkg.add("JuMP", "CSV", "DataFrames", "PlotlyJS", "Dates", "XLSX", "FileIO", "Base", "Random", "Statistics", "Gurobi", "PyCall")
    # Pkg.build("PyCall")
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
    include("User_Input.jl")
    include("Input_Parameters.jl")
    include("Data_Conversion.jl")
    include("PassiveModel_5.6.jl")
    include("SolarGen V1.0.jl")
    include("ElectricLoad V1.0.jl")
    include("RuleBased.jl")

    ENV["PYTHON"] = Python_Environment
    ENV["GRB_LICENSE_FILE"] = Gurobi_License

    pushfirst!(PyVector(pyimport("sys")."path"), "C:\\Users\\Fred\\Desktop\\PyFMI")

    ############ Load FMU ############
    begin
        pyf = pyimport("pyfmi")
        FMU = pyimport(FMU_Simulation_File)
  
        # Load the FMU model
        model = pyf.load_fmu(fmu_Name)
    end
end

############ Program Preparations ############
begin
    # Update automatically the date when this program is run.
    today_date = today()

    # Please update information of this program to automatically update the code name.
    version = 1.0

    # Create folder to later save data and plots
    begin
        # Define the path for the new folder
        folder_path = joinpath(Main_Folder_Path, "$code_name")

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
    # Check if the file exists before trying to read it
    if isfile(weather_file_name)
        # weather = DataFrame(CSV.File(file_name, header=3))
        weather = DataFrame(CSV.File(weather_file_name))
    else
        println("File $weather_file_name not found.")
    end
end

############ Declare Parameters ############
# Time Horizon Parameters
begin
    stepsize = 60*(1/f_run) # [min]
    # base_stepsize = 5 # [min]
    # base_f_run = 60/base_stepsize
    # TimeEnd = 8760 * f_run;
    #=
    steps = Int[]
    for i in eachindex(stepsizes)
        append!(steps, fill(stepsizes[i], stepnums[i]))
    end
    =#
end    

df = convert_weather(weather, "epw")
Q_total, Q_cd, Q_cv, Q_rh, Q_rc, TC = passive_model(calibration_file_path, df, T_indoor_constant)
pv_cf = Generate_PV(df)
e_load = generate_schedules("complex")

df = combine_dfs(df, Q_total, pv_cf, e_load)

input_df = prepare_input_rule_based(df)

############ Initialize Space ############
begin
    # Initialize M_States to start the FMU simulation and Julia optimization
    #=
    These are the current initial conditions:
        PCM Hot Initial Temperature = 48 [°C] (half full)
        PCM Cold Initial Temperature = 11 [°C] (half full)
        Initial Indoor Temperature = 25.5 [°C]
        PCM Hot Capacity = PCM_H_Size = 9 [kWh]
        PCM Cold Capacity = PCM_C_Size = 9 [kWh]
        Zone Temperature Cooling Setpoint = 28 [°C]
        Zone Temperature Heating Setpoint = 23 [°C]
        Zone Temperature Delta Setpoint = 2 [°C]
    =#

    global current_m_states = DataFrame(PCM_Hot_Temp = [273.15 + Inital_PCM_H_Temp], PCM_Cold_Temp = [273.15 + Inital_PCM_C_Temp], Indoor_Temp = [273.15 + Initial_Ti], 
    PCM_Hot_Cap = [PCM_H_Size], PCM_Cold_Cap = [PCM_C_Size],  
    T_Cooling = [273.15 + zone_temp_cooling_setpoint], T_Heating = [273.1 + zone_temp_heating_setpoint], 
    T_Delta = [zone_temp_setpoint_delta]);
    
    # This is the same as current_m_states at i=1, but it works with the python input for the FMU.
    initialstates = [];
    for i = 1:size(current_m_states)[2]
        push!(initialstates, current_m_states[1, i])
    end

    # Initialize J_States: [PV_generation, E_0, InStorageBattery]
    #=
    These are the current initial conditions:
        PV Generation at last timestep = 0 [kWh]
        Actual House load consumption at last timestep = 0 [kWh]
        Battery Initial SOC = 0.5
    =#

    initial_j_states = DataFrame(
        "PV Generation (kW)" => [Initial_PV_Gen],
        "Lighting Plug Load (kW)" => [Initial_P_0],
        "Remain Battery Energy (kWh)" => [Initial_B_SOC * BatterySize],
        "Curtailment (kW)" => [Initial_Curtailment],
        "Loss of Load (kW)"  => [Initial_Loss_of_Load]
    )
    J_States_list = [initial_j_states]

    Actions = zeros(NumRun);
    Costs_e = zeros(NumRun);
    
    Julia_Scenario = zeros(NumRun);

    fmu_temp = zeros(NumRun);
    fmu_PCM_C_SOC = zeros(NumRun);
    fmu_PCM_H_SOC = zeros(NumRun);
    fmu_indoor_temp = zeros(NumRun);
    fmu_Pump1 = zeros(NumRun);
    fmu_Pump2 = zeros(NumRun);
    fmu_HP_Mode = zeros(NumRun);
    fmu_HP_ON_OFF = zeros(NumRun);
    fmu_HP_power = zeros(NumRun);
    fmu_Useful_Thermal= zeros(NumRun);
    fmu_fan_power = zeros(NumRun);
    fmu_PCM_C_Temp = zeros(NumRun);
    fmu_PCM_H_Temp = zeros(NumRun);
    fmu_P1_Mass_Flow = zeros(NumRun);
    fmu_P2_Mass_Flow = zeros(NumRun);
    fmu_HT_PCM_C = zeros(NumRun);
    fmu_HT_PCM_H = zeros(NumRun);
    fmu_FCU_Heat = zeros(NumRun);
    fmu_system_Mode = zeros(NumRun);
    fmu_FCU_W_Enter_Temp = zeros(NumRun);
    fmu_FCU_W_Leave_Temp = zeros(NumRun);
    fmu_Fan_Coil_Mass_Flow = zeros(NumRun);
    fmu_HP_W_Supply_Temp = zeros(NumRun);
    fmu_HP_W_Return_Temp = zeros(NumRun);
    fmu_HP_Compressor_Speed = zeros(NumRun);
    fmu_FCU_Mode = zeros(NumRun);
    fmu_Mode_from_MPC = zeros(NumRun);
    fmu_PCM_H_empty = zeros(NumRun);
    fmu_PCM_C_empty = zeros(NumRun);
    fmu_HP_W_Supply_Temp_SP = zeros(NumRun);

    runtime = zeros(NumRun);
    # Heat_pump_modes_j = zeros(NumRun);
    # HP_H_Usage_j = zeros(NumRun);
end

############ Program Execution ############
runtime1 = @elapsed begin
    starttime = weather.DateTime[TimeStart];
    
    initial_timepoints = [0, 60*stepsize]
    initial_schedule = [0, 0]; # Do nothing while the model is solving the first iteration

    res = FMU.fmu(1, initialstates, initial_schedule, initial_timepoints, model)
    initial_m_states = pandas_to_julia(res)
    
    initial_m_states = convert_temperatures(initial_m_states) # temperatures now in C

    fmu_temp[1] = initial_m_states[end, "Ambient Temp"];
    fmu_PCM_C_SOC[1], fmu_PCM_H_SOC[1] = PCM_available_energy(initial_m_states[end, "PCM_Cold_Temp"], initial_m_states[end, "PCM_Hot_Temp"]);
    fmu_indoor_temp[1] = initial_m_states[end, "Indoor_Temp"];
    fmu_HP_Mode[1] = initial_m_states[end, "HP Mode"];
    fmu_HP_ON_OFF[1] = initial_m_states[end, "HP ON OFF"];
    fmu_PCM_H_empty[1] = initial_m_states[end, "hot pcm not fully discharged"];
    fmu_PCM_C_empty[1] = initial_m_states[end, "cold pcm not fully discharged"];

    fmu_Mode_from_MPC[1] = initial_m_states[end, "FMU mode from MPC"];
    fmu_FCU_Mode[1] = initial_m_states[end, "FCU Mode"];
    fmu_system_Mode[1] = initial_m_states[end, "Exact System Mode"];
    fmu_PCM_H_Temp[1] = initial_m_states[end, "PCM_Hot_Temp"];
    fmu_PCM_C_Temp[1] = initial_m_states[end, "PCM_Cold_Temp"];
    fmu_FCU_W_Enter_Temp[1] = initial_m_states[end, "FCU Entering Water Temperature"];
    fmu_FCU_W_Leave_Temp[1] = initial_m_states[end, "FCU Leaving Water Temperature"];
    fmu_HP_W_Supply_Temp[1] = initial_m_states[end, "Heat Pump Supply Water Temperature"];
    fmu_HP_W_Return_Temp[1] = initial_m_states[end, "Heat Pump Return Water Temperature"];
    fmu_HP_W_Supply_Temp_SP[1] = initial_m_states[end, "Heat Pump Supply Water Temperature Setpoint"];

    # Average Rate 
    fmu_Fan_Coil_Mass_Flow[1] = find_average(initial_m_states, "Fan Coil Fan Mass Flow"); #[kg/s]
    fmu_HP_Compressor_Speed[1] = find_average(initial_m_states, "Heat Pump Compressor Speed");
    fmu_P1_Mass_Flow[1] = find_average(initial_m_states, "Pump 1 Mass Flow"); #[kg/s]
    fmu_P2_Mass_Flow[1] = find_average(initial_m_states, "Pump 2 Mass Flow"); #[kg/s]

    # Average Power
    fmu_fan_power[1] = p2e(initial_m_states)[1]*f_run; #[kW]
    fmu_FCU_Heat[1] = p2e(initial_m_states)[2]*f_run; #[kW]
    fmu_HP_power[1]= p2e(initial_m_states)[3]*f_run; #[kW]
    fmu_HT_PCM_C[1] = p2e(initial_m_states)[4]*f_run; #[kW]
    fmu_HT_PCM_H[1] = p2e(initial_m_states)[5]*f_run; #[kW]
    fmu_Pump1[1] = p2e(initial_m_states)[6]*f_run; #[kW]
    fmu_Pump2[1] = p2e(initial_m_states)[7]*f_run; #[kW]
    fmu_Useful_Thermal[1]= p2e(initial_m_states)[8]*f_run; #[kW]

    # Updates the Julia and Modelica states.
    global current_m_states = initial_m_states
    global current_j_states = initial_j_states

    M_States_list = [initial_m_states]
end
runtime[1]= runtime1; #[seconds]

global previous_decision = 0

for i = TimeStart+1:1:NumRun
    current_runtime = @elapsed begin
        println()
        println("Iteration: ", i)
        println()

        # Define timepoints needed for FMU simulation
        timepoints = [Dates.value(weather.DateTime[i]-starttime)/1000, Dates.value(weather.DateTime[i+1]-starttime)/1000] # [s]

        lol_e, new_j_states, input_schedule, scenario = RuleBased(i-1, input_df[i:i+f_run*120,:], current_m_states, current_j_states, previous_decision);
        
        global previous_decision = input_schedule[1]

        println("Successful!")
        println()
        # println("Input schedule: ", input_schedule)
        # println()
        println("timepoints: ", timepoints)
        println()
        println("Calling FMU")

        # Store actions
        Actions[i] = input_schedule[1]

        # Add the lost of load from the current iteration to total lost of load (cost)
        Costs_e[i] = lol_e; # [kW]

        # Call the FMU function with the input (current iteration, operational commands, timepoints, FMU model)
        res = FMU.fmu(i, initialstates, input_schedule, timepoints, model)
        
        # This is the dataframe output from the current FMU iteration with power consumption and SOC updates etc. (Modelica states)
        new_m_states = pandas_to_julia(res)

        new_m_states = convert_temperatures(new_m_states) # temperatures now in C

        Julia_Scenario[i] = scenario # the command for the current iteration 
        
        # These are the information at the end of the current iteration, which become input for the start of the next iteration
        fmu_temp[i] = new_m_states[end, "Ambient Temp"];
        fmu_PCM_C_SOC[i], fmu_PCM_H_SOC[i] = PCM_available_energy(new_m_states[end, "PCM_Cold_Temp"], new_m_states[end, "PCM_Hot_Temp"]);
        fmu_indoor_temp[i] = new_m_states[end, "Indoor_Temp"];
        fmu_HP_Mode[i] = new_m_states[end, "HP Mode"];
        fmu_HP_ON_OFF[i] = new_m_states[end, "HP ON OFF"];

        fmu_PCM_H_empty[i] = new_m_states[end, "hot pcm not fully discharged"];
        fmu_PCM_C_empty[i] = new_m_states[end, "cold pcm not fully discharged"];
        fmu_Mode_from_MPC[i] = new_m_states[end, "FMU mode from MPC"];
        fmu_FCU_Mode[i] = new_m_states[end, "FCU Mode"];
        fmu_system_Mode[i] = new_m_states[end, "Exact System Mode"];
        fmu_PCM_H_Temp[i] = new_m_states[end, "PCM_Hot_Temp"];
        fmu_PCM_C_Temp[i] = new_m_states[end, "PCM_Cold_Temp"];
        fmu_FCU_W_Enter_Temp[i] = new_m_states[end, "FCU Entering Water Temperature"];
        fmu_FCU_W_Leave_Temp[i] = new_m_states[end, "FCU Leaving Water Temperature"];
        fmu_HP_W_Supply_Temp[i] = new_m_states[end, "Heat Pump Supply Water Temperature"];
        fmu_HP_W_Return_Temp[i] = new_m_states[end, "Heat Pump Return Water Temperature"];
        fmu_HP_W_Supply_Temp_SP[i] = new_m_states[end, "Heat Pump Supply Water Temperature Setpoint"];

        # Average Rate 
        fmu_Fan_Coil_Mass_Flow[i] = find_average(new_m_states, "Fan Coil Fan Mass Flow"); #[kg/s]
        fmu_HP_Compressor_Speed[i] = find_average(new_m_states, "Heat Pump Compressor Speed");
        fmu_P1_Mass_Flow[i] = find_average(new_m_states, "Pump 1 Mass Flow"); #[kg/s]
        fmu_P2_Mass_Flow[i] = find_average(new_m_states, "Pump 2 Mass Flow"); #[kg/s]

        # Average Power
        fmu_fan_power[i] = p2e(new_m_states)[1]*f_run; #[kW]
        fmu_FCU_Heat[i] = p2e(new_m_states)[2]*f_run; #[kW]
        fmu_HP_power[i]= p2e(new_m_states)[3]*f_run; #[kW]
        fmu_HT_PCM_C[i] = p2e(new_m_states)[4]*f_run; #[kW]
        fmu_HT_PCM_H[i] = p2e(new_m_states)[5]*f_run; #[kW]
        fmu_Pump1[i] = p2e(new_m_states)[6]*f_run; #[kW]
        fmu_Pump2[i] = p2e(new_m_states)[7]*f_run; #[kW]
        fmu_Useful_Thermal[i]= p2e(new_m_states)[8]*f_run; #[kW]

        # Save the Julia and Modelica states.
        push!(M_States_list, new_m_states);
        push!(J_States_list, new_j_states);

        # Updates the Julia and Modelica states.
        global current_m_states = new_m_states;
        global current_j_states = new_j_states;
    end
    runtime[i] = current_runtime # [s]
end


pop!(Actions)
insert!(Actions, 1, 0)
push!(Actions, 0)

insert!(Julia_Scenario, 1, 0)

push!(Costs_e, 0)
push!(Costs_q, 0)
push!(runtime, 0)

# These M_States are at the start of the 00:00, therefore they were inserted to index 1

fmu_PCM_C_SOC_1, fmu_PCM_H_SOC_1 = PCM_available_energy(initial_m_states[1, "PCM_Cold_Temp"], initial_m_states[1, "PCM_Hot_Temp"]);
insert!(fmu_PCM_H_SOC, 1, fmu_PCM_H_SOC_1)
insert!(fmu_PCM_C_SOC, 1, fmu_PCM_C_SOC_1)
insert!(fmu_temp, 1, initial_m_states[1, "Ambient Temp"])
insert!(fmu_indoor_temp, 1, initial_m_states[1, "Indoor_Temp"])
insert!(fmu_HP_Mode, 1, initial_m_states[1, "HP Mode"])
insert!(fmu_HP_ON_OFF, 1, initial_m_states[1, "HP ON OFF"])
insert!(fmu_FCU_Mode, 1, initial_m_states[1, "FCU Mode"])
insert!(fmu_Mode_from_MPC, 1, initial_m_states[1, "FMU mode from MPC"])
insert!(fmu_system_Mode, 1, initial_m_states[1, "Exact System Mode"])
insert!(fmu_PCM_H_Temp, 1, initial_m_states[1, "PCM_Hot_Temp"])
insert!(fmu_PCM_C_Temp, 1, initial_m_states[1, "PCM_Cold_Temp"])
insert!(fmu_FCU_W_Enter_Temp, 1, initial_m_states[1, "FCU Entering Water Temperature"])
insert!(fmu_FCU_W_Leave_Temp, 1, initial_m_states[1, "FCU Leaving Water Temperature"])
insert!(fmu_HP_W_Supply_Temp, 1, initial_m_states[1, "Heat Pump Supply Water Temperature"])
insert!(fmu_HP_W_Return_Temp, 1, initial_m_states[1, "Heat Pump Return Water Temperature"])
insert!(fmu_HP_W_Supply_Temp_SP, 1, initial_m_states[1, "Heat Pump Supply Water Temperature Setpoint"])
insert!(fmu_PCM_H_empty, 1, initial_m_states[1, "hot pcm not fully discharged"])
insert!(fmu_PCM_C_empty, 1, initial_m_states[1, "cold pcm not fully discharged"])

# These M_states are arbitrary, they are filler values. They are called flow states.
# Let's say we do 2 iterations in total:
# First iteration: start at 00:00, end at 00:05
# Second iteration: start at 00:05, end at 00:10
# We still need to add an extra row that starts at 00:10 to record the end M_states at 00:10, but since there is no process between 00:10 to 00:15, the flow states do not exist.

push!(fmu_Fan_Coil_Mass_Flow, 0)
push!(fmu_HP_Compressor_Speed, 0)
push!(fmu_P1_Mass_Flow, 0)
push!(fmu_P2_Mass_Flow, 0)
push!(fmu_fan_power, 0)
push!(fmu_FCU_Heat, 0)
push!(fmu_HP_power, 0)
push!(fmu_HT_PCM_C, 0)
push!(fmu_HT_PCM_H, 0)
push!(fmu_Pump1, 0)
push!(fmu_Pump2, 0)
push!(fmu_Useful_Thermal, 0)

# Hence the reason why we extend the Datetime array by an extra index.
Datetime_array = input_df[1:NumRun+1, "DateTime"]

# The PV, Lighting & Plug loads are the same as flow states. They are filler values in the last row.

# Battery SOC should be solved at the extra timestep using the initial battery SOC from the very last timestep and PV and total load during the very last timestep.

# Curtailment and Loss of Load are flow states from the previous iteration. They should be solved at the extra timestep to get the last timestep, and then include filler values.

final_j_states = DataFrame(
    "PV Generation (kW)" => [0],
    "Lighting Plug Load (kW)" => [0],
    "Remain Battery Energy (kWh)" => [0],
    "Curtailment (kW)" => [0],
    "Loss of Load (kW)"  => [0]
)

push!(J_States_list, final_j_states)

M_df_combined = reduce(vcat, M_States_list)
J_df_combined = reduce(vcat, J_States_list)

# EXPORT INTO A CSV FILE (Time series)

############ Save Results to Excel ############
begin
    # Create Excel files
    begin
        # Record the exact time when the result is produced
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        # Create file names and sheet names
        xlsx_file = "$code_name($timestamp).xlsx"
        sheetnames = ["Discrete Data", "Continuous Data"]
    end
    
    # Store Result Data
    begin
        column_names = ["DateTime", "Run Time (s)", "Ti (C)", "Ta (C)", "MPC Command", "FMU Mode Received from MPC", "Julia Scenario", "FMU System Mode", "FMU FCU Mode", "FMU HP Mode", "FMU HP on/off",
        "PCM H Remaining Energy (kWh)", "PCM C Remaining Energy (kWh)", "PCM H Temp (C)", "PCM C Temp (C)", "PCM H Not Fully Discharged", "PCM C Not Fully Discharged", 
        "PCM H Ave Heat Transfer (kW)", "PCM C Ave Heat Transfer (kW)", "Useful Ave Thermal Delivered (kW)",
        "FCU Entering Water Temp", "FCU Leaving Water Temp", "Heat Pump Supply Water Temp", "Heat Pump Return Water Temp", "Heat Pump Supply Water Temp Setpoint",
        "FCU Ave Heat Power (kW)", "Fan Ave Power (kW)", "HP Ave Power (kW)", 
        "Pump 1 Ave Power (kW)", "Pump 2 Ave Power (kW)", "Heat Pump Compressor Speed", 
        "Fan Coil Fan Ave Mass Flow (kg/s)", "Pump 1 Ave Mass Flow (kg/s)", "Pump 2 Ave Mass Flow (kg/s)", 
        "Loss of Load (Electrical) (kW)", "Loss of Load (Thermal) (kW)"]

        # Collect lists into a tuple or an array
        data_lists = [Datetime_array, runtime, fmu_indoor_temp, fmu_temp, Actions, fmu_Mode_from_MPC, Julia_Scenario, fmu_system_Mode, fmu_FCU_Mode, fmu_HP_Mode, fmu_HP_ON_OFF,
        fmu_PCM_H_SOC, fmu_PCM_C_SOC, fmu_PCM_H_Temp, fmu_PCM_C_Temp, fmu_PCM_H_empty, fmu_PCM_C_empty,
        fmu_HT_PCM_H, fmu_HT_PCM_C, fmu_Useful_Thermal, 
        fmu_FCU_W_Enter_Temp, fmu_FCU_W_Leave_Temp, fmu_HP_W_Supply_Temp, fmu_HP_W_Return_Temp, fmu_HP_W_Supply_Temp_SP,
        fmu_FCU_Heat, fmu_fan_power, fmu_HP_power,
        fmu_Pump1, fmu_Pump2, fmu_HP_Compressor_Speed, 
        fmu_Fan_Coil_Mass_Flow, fmu_P1_Mass_Flow, fmu_P2_Mass_Flow,
        Costs_e, Costs_q]
        
        Discrete_Data = DataFrame(data_lists, column_names)
        Discrete_Data = hcat(Discrete_Data, J_df_combined)
    end
    # Define the path for the Excel file
    excel_file_path = joinpath(folder_path, xlsx_file)

    # Write the DataFrames to different sheets in the same Excel file
    XLSX.openxlsx(excel_file_path, mode="w") do xf

        # Add and write to the first sheet
        sheet1 = XLSX.addsheet!(xf, sheetnames[1])
        XLSX.writetable!(sheet1, DataFrames.eachcol(Discrete_Data), DataFrames.names(Discrete_Data))

        # Add and write to the second sheet
        sheet2 = XLSX.addsheet!(xf, sheetnames[2])
        XLSX.writetable!(sheet2, DataFrames.eachcol(M_df_combined), DataFrames.names(M_df_combined))
    end
end
