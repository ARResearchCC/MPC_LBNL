# Model Predictive Control: LBNL Simulation. Version 5.0
# Complete Date: 10/27/2024
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
    include("MPCAlgorithms_LBNL_New.jl")

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
    version = 5.0

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
    base_stepsize = 5 # [min]
    base_f_run = 60/base_stepsize
    # TimeEnd = 8760 * f_run;

    steps = Int[]
    for i in eachindex(stepsizes)
        append!(steps, fill(stepsizes[i], stepnums[i]))
    end
end    

############ Initialize Space ############
begin
    # Initialize M_States to start the FMU simulation and Julia optimization
    #=
    These are the current initial conditions:
        PCM Hot Initial Temperature = 48 [°C]
        PCM Cold Initial Temperature = 11 [°C]
        Initial Indoor Temperature = 22 [°C]
        PCM Hot Capacity = PCM_H_Size = 11 [kWh]
        PCM Cold Capacity = PCM_C_Size = 7 [kWh]
        Zone Temperature Cooling Setpoint = 22 [°C]
        Zone Temperature Heating Setpoint = 20 [°C]
        Zone Temperature Delta Setpoint = 0.5 [°C]
        PCM Cold Initial SOC = 0.5
        PCM Hot Initial SOC = 0.5
    =#

    global current_m_states = DataFrame(PCM_Hot_Temp = [273.15 + Inital_PCM_H_Temp], PCM_Cold_Temp = [273.15 + Inital_PCM_C_Temp], Indoor_Temp = [273.15 + Initial_Ti], 
    PCM_Hot_Cap = [PCM_H_Size], PCM_Cold_Cap = [PCM_C_Size],  
    T_Cooling = [273.15 + zone_temp_cooling_setpoint], T_Heating = [273.1 + zone_temp_heating_setpoint], 
    T_Delta = [zone_temp_setpoint_delta], PCM_Cold_SOC = [Intial_PCM_C_SOC], PCM_Hot_SOC = [Intial_PCM_H_SOC]);
    
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
        "PV Generation (kWh)" => [Initial_PV_Gen],
        "Lighting Plug Load (kWh)" => [Initial_P_0],
        "Remain Battery Energy (kWh)" => [Initial_B_SOC * BatterySize],
        "Curtailment (kWh)" => [Initial_Curtailment]
    )
    J_States_list = [initial_j_states]

    Solve_failure = zeros(NumRun);
    Actions = zeros(NumRun);
    Costs = zeros(NumRun);
    AccumulatedCosts = zeros(NumRun);
    global TotalCost = 0;

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

############ Data Preparations ############

df = convert_weather(weather, "epw")
Q_total, Q_cd, Q_cv, Q_rh, Q_rc, TC = passive_model(calibration_file_path, df, T_indoor_constant)
pv_cf = Generate_PV(df)
e_load = generate_schedules("complex")

df = combine_dfs(df, Q_total, pv_cf, e_load)

############ Program Execution ############
runtime1 = @elapsed begin
    starttime = weather.DateTime[TimeStart];
    
    initial_timepoints = [0, 60*stepsize]
    initial_schedule = [0, 0]; # Do nothing while the model is solving the first iteration
    # Actions[1] = 0;

    res = FMU.fmu(1, initialstates, initial_schedule, initial_timepoints, model)
    initial_m_states = pandas_to_julia(res)
    
    fmu_temp[1] = initial_m_states[end, "Ambient Temp"];
    fmu_PCM_H_SOC[1] = initial_m_states[end, "PCM_Hot_SOC"];
    fmu_PCM_C_SOC[1] = initial_m_states[end, "PCM_Cold_SOC"];
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

for i = TimeStart+1:1:NumRun
    current_runtime = @elapsed begin
        println()
        println("Iteration: ", i)
        println()
        #=
        println("The current M states are:")
        println(current_m_states)
        =#

        # Define timepoints needed for FMU simulation
        timepoints = [Dates.value(weather.DateTime[i]-starttime)/1000, Dates.value(weather.DateTime[i+1]-starttime)/1000] # [s]
        # println(timepoints)

        # input_data = resample_mpc(df, i, steps, stepsize)
        input_data = aggregate_energy(df, i, steps, base_stepsize)
        
        # println("The input data from the passive model is:")
        # println(input_data)
        # Call the MPC function
        currentcost, new_j_states, input_schedule, didnt_solve = OptimizeV4_1(i-1, TC, steps, base_stepsize, input_data, current_m_states, current_j_states);
        
        println("Optimization Successful!")
        println()
        println("Input schedule: ", input_schedule)
        println()
        println("timepoints: ", timepoints)
        println()
        println("Calling FMU")

        # Is the current iteration of optimization infeasible?
        Solve_failure[i] = didnt_solve
        # Store actions
        Actions[i] = input_schedule[1]
        # Add the lost of load (cost) from the current iteration to total lost of load (cost)
        Costs[i] = currentcost;
        AccumulatedCosts[i] = sum(Costs[n] for n = 1:i)
        global TotalCost += currentcost;

        # Call the FMU function with the input (current iteration, operational commands, timepoints, FMU model)
        res = FMU.fmu(i, initialstates, input_schedule, timepoints, model)
        # println()
        # println("The Modelica Feedbacks are:")
        # println(res)
        
        # This is the dataframe output from the current FMU iteration with power consumption and SOC updates etc. (Modelica states)
        new_m_states = pandas_to_julia(res)

        fmu_temp[i] = new_m_states[end, "Ambient Temp"];
        fmu_PCM_H_SOC[i] = new_m_states[end, "PCM_Hot_SOC"];
        fmu_PCM_C_SOC[i] = new_m_states[end, "PCM_Cold_SOC"];
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

final_j_states = DataFrame(
    "PV Generation (kWh)" => [0],
    "Lighting Plug Load (kWh)" => [0],
    "Remain Battery Energy (kWh)" => [0],
    "Curtailment (kWh)" => [0]
)

push!(J_States_list, final_j_states)
pop!(Actions)
insert!(Actions, 1, 0)
push!(Actions, 0)
push!(Costs, 0)
push!(AccumulatedCosts, AccumulatedCosts[end])
push!(Solve_failure, 0)
push!(runtime, 0)

insert!(fmu_PCM_H_empty, 1, initial_m_states[1, "hot pcm not fully discharged"])
insert!(fmu_PCM_C_empty, 1, initial_m_states[1, "cold pcm not fully discharged"])
insert!(fmu_temp, 1, initial_m_states[1, "Ambient Temp"])
insert!(fmu_PCM_H_SOC, 1, initial_m_states[1, "PCM_Hot_SOC"])
insert!(fmu_PCM_C_SOC, 1, initial_m_states[1, "PCM_Cold_SOC"])
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

# Select the desired rows
resample_ratio = div(stepsize, base_stepsize)
selected_rows = Int.(collect(1:resample_ratio:size(df, 1)))
resampled_df = df[selected_rows, :]
Datetime_array = resampled_df[1:NumRun+1, "DateTime"]
begin
#=
    #=
    # Plot Accumulated loss of load
    trace1 = scatter(
        x = df[1:NumRun, "DateTime"],  
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
        x = df[1:NumRun, "DateTime"],  
        y=fmu_temp.-fahrenheit_to_kelvin(df[1:NumRun, "Ta"]),
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
        x = df[1:NumRun, "DateTime"],  
        y=fmu_indoor_temp.-273.15,
        name="Indoor Temp"
    )

    trace9 = scatter(
        x = df[1:NumRun, "DateTime"],  
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
        x = df[1:NumRun, "DateTime"],  
        y=AllStates[:, 1]/BatterySize,
        name="Battery SOC"
    )
    trace5 = scatter(
        x = df[1:NumRun, "DateTime"],  
        y=AllStates[:, 2],
        name="PCM H SOC"
    )
    trace6 = scatter(
        x = df[1:NumRun, "DateTime"],  
        y=AllStates[:, 3],
        name="PCM C SOC"
    )

    p_soc = plot([trace4, trace5, trace6], Layout(title="SOC vs. Time", xaxis_title="Time Steps", yaxis_title="SOC percentage"))
    display(p_soc)

    trace10 = scatter(
        x = df[1:NumRun, "DateTime"],  
        y=fmu_fan,
        name="FMU Fan Energy Consumption"
    )

    trace11 = scatter(
        x = df[1:NumRun, "DateTime"],  
        y=fmu_FCU,
        name="FMU Fan Coil Unit Heat Delivered"
    )

    trace13 = scatter(
        x = df[1:NumRun, "DateTime"],  
        y=fmu_HP_power,
        name="fmu HP Energy Consumption"
    )

    trace14 = scatter(
        x = df[1:NumRun, "DateTime"],  
        y=fmu_HT_PCM_C,
        name="FMU Heat Transfer PCM C"
    )

    trace15 = scatter(
        x = df[1:NumRun, "DateTime"],  
        y=fmu_HT_PCM_H,
        name="FMU Heat Transfer PCM H"
    )

    trace16 = scatter(
        x = df[1:NumRun, "DateTime"],  
        y=fmu_Pump1,
        name="FMU Pump 1"
    )

    trace17 = scatter(
        x = df[1:NumRun, "DateTime"],  
        y=fmu_Pump2,
        name="FMU Pump 2"
    )

    trace18 = scatter(
        x = df[1:NumRun, "DateTime"],  
        y=fmu_Useful_Thermal,
        name="FMU Useful Thermal Delivered"
    )

    trace19 = scatter(
        x = df[1:NumRun, "DateTime"],  
        y=AllStates[:, 2].*(PCM_H_Size),
        name="PCM H Remaining Energy"
    )
    trace20 = scatter(
        x = df[1:NumRun, "DateTime"],  
        y=AllStates[:, 3].*(PCM_C_Size),
        name="PCM C Remaining Energy"
    )

    trace21 = scatter(
        x = df[1:NumRun, "DateTime"],  
        y=Curtailments,
        name="PV Curtailed Energy"
    )

    p_energy = plot([trace10, trace16, trace17], Layout(title="Device Energy Consumption every 5 minutes", xaxis_title="Time Steps", yaxis_title="Energy Consumption (kWh)"))
    display(p_energy)

    p_thermal = plot([trace18, trace11, trace13], Layout(title="Device Thermal Energy every 5 minutes", xaxis_title="Time Steps", yaxis_title="Energy (kWh)"))
    display(p_thermal)

    p_pcm = plot([trace14, trace15, trace19, trace20], Layout(title="PCM Plots", xaxis_title="Time Steps", yaxis_title="Energy (kWh)"))
    display(p_pcm)

    p_curtailed = plot([trace21], Layout(title="Curtailment Plots", xaxis_title="Time Steps", yaxis_title="Energy (kWh)"))
    display(p_curtailed)
=#    
end

# Runtime per iteration (must be less than the MPC run frequency)

println("The total loss of load is $TotalCost kWh")

Total_Inf = sum(Solve_failure)
Inf_pct = 100*Total_Inf/(NumRun-1)

println("There are a total of $Total_Inf infeasible iterations out of $(NumRun-1) iterations, which is $Inf_pct % of the total iterations.")

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
        xlsx_file = "MPC($version)($timestamp).xlsx"
        sheetnames = ["Discrete Data", "Continuous Data"]
    end
    
    # Store Result Data
    begin
        column_names = ["DateTime", "Run Time (s)", "Julia Solve Failure", "Ti (C)", "Ta (C)", "MPC Command", "FMU Mode Received from MPC", "FMU System Mode", "FMU FCU Mode", "FMU HP Mode", "FMU HP on/off",
        "PCM H SOC", "PCM C SOC", "PCM H Remaining Energy (kWh)", "PCM C Remaining Energy (kWh)", "PCM H Temp (C)", "PCM C Temp (C)", "PCM H Not Fully Discharged", "PCM C Not Fully Discharged", 
        "PCM H Ave Heat Transfer (kW)", "PCM C Ave Heat Transfer (kW)", "Useful Ave Thermal Delivered (kW)",
        "FCU Entering Water Temp", "FCU Leaving Water Temp", "Heat Pump Supply Water Temp", "Heat Pump Return Water Temp", "Heat Pump Supply Water Temp Setpoint",
        "FCU Ave Heat Power (kW)", "Fan Ave Power (kW)", "HP Ave Power (kW)", 
        "Pump 1 Ave Power (kW)", "Pump 2 Ave Power (kW)", "Heat Pump Compressor Speed", 
        "Fan Coil Fan Ave Mass Flow (kg/s)", "Pump 1 Ave Mass Flow (kg/s)", "Pump 2 Ave Mass Flow (kg/s)", 
        "Loss of Load (kWh)", "Total Loss of Load (kWh)"]

        # Collect lists into a tuple or an array
        data_lists = [Datetime_array, runtime, Solve_failure, fmu_indoor_temp.-273.15, fmu_temp.-273.15, Actions, fmu_Mode_from_MPC, fmu_system_Mode, fmu_FCU_Mode, fmu_HP_Mode, fmu_HP_ON_OFF,
        fmu_PCM_H_SOC, fmu_PCM_C_SOC, fmu_PCM_H_SOC.*(PCM_H_Size), fmu_PCM_C_SOC.*(PCM_C_Size), fmu_PCM_H_Temp.-273.15, fmu_PCM_C_Temp.-273.15, fmu_PCM_H_empty, fmu_PCM_C_empty,
        fmu_HT_PCM_H, fmu_HT_PCM_C, fmu_Useful_Thermal, 
        fmu_FCU_W_Enter_Temp.-273.15, fmu_FCU_W_Leave_Temp.-273.15, fmu_HP_W_Supply_Temp.-273.15, fmu_HP_W_Return_Temp.-273.15, fmu_HP_W_Supply_Temp_SP.-273.15,
        fmu_FCU_Heat, fmu_fan_power, fmu_HP_power,
        fmu_Pump1, fmu_Pump2, fmu_HP_Compressor_Speed, 
        fmu_Fan_Coil_Mass_Flow, fmu_P1_Mass_Flow, fmu_P2_Mass_Flow,
        Costs, AccumulatedCosts]
        
        Discrete_Data = DataFrame(data_lists, column_names)
        Discrete_Data = hcat(Discrete_Data, J_df_combined)
    end
    # Define the path for the Excel file
    excel_file_path = joinpath(folder_path, xlsx_file)

    # Write the DataFrames to different sheets in the same Excel file
    XLSX.openxlsx(excel_file_path, mode="w") do xf
        sheet1 = XLSX.addsheet!(xf, sheetnames[1])
        XLSX.writetable!(sheet1, DataFrames.eachcol(Discrete_Data), DataFrames.names(Discrete_Data))

        sheet2 = XLSX.addsheet!(xf, sheetnames[2])
        XLSX.writetable!(sheet2, DataFrames.eachcol(M_df_combined), DataFrames.names(M_df_combined))
    end
end
