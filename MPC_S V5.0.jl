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
    # Pkg.update("Gurobi")
    # ENV["PYTHON"] = "C:\\Users\\Fred\\anaconda3\\python.exe"
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
    include("Input_Parameters.jl")
    include("Data_Conversion.jl")
    include("PassiveModel_5.6.jl")
    include("SolarGen V1.0.jl")
    include("ElectricLoad V1.0.jl")
    include("MPCAlgorithms_LBNL_New.jl")


    ENV["GRB_LICENSE_FILE"] = "C:\\Users\\Fred\\gurobi.lic" # Fred's license

    pushfirst!(PyVector(pyimport("sys")."path"), "C:\\Users\\Fred\\Desktop\\PyFMI")

    ############ Load FMU ############
    begin
        pyf = pyimport("pyfmi")
        FMU = pyimport("FMU_Simulation_4")
        # Define the path to the FMU file
        # fmuName = "Stanford_Hybrid_System.fmu"
  
        # Load the FMU model
        model = pyf.load_fmu("C:\\Users\\Fred\\Desktop\\PyFMI\\Stanford_Hybrid_System.fmu")
    end
end

############ Program Preparations ############
begin
    # Update automatically the date when this program is run.
    today_date = today()

    # Please update information of this program to automatically update the code name.
    version = 5.0

    code_name = "MPC_S_V$version._$today_date"

    # Create folder to later save data and plots
    begin
        # Define the path for the new folder
        folder_path = "C:\\Users\\Fred\\Desktop\\PyFMI\\$code_name"

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

    calibration_filepath = "C:\\Users\\Fred\\Desktop\\PyFMI\\Calibration_Model_Input.xlsx"
end

############ Declare Parameters ############
# Time Horizon Parameters
begin
    TimeStart = 1;
    f_run = 12 # run frequency [times/hour]
    stepsize = 60*(1/f_run) # [min]
    TimeEnd = 8760 * f_run;
    # Opt_Horizon = 6 * f_run # [intervals]
    # NumRun = (TimeEnd-TimeStart+1) - Opt_Horizon + 1; # This is the max number of runs allowed for a full year of data
    NumRun = 24*12

    # Define the sampling intervals
    stepsizes = [5, 30, 60, 120] # minutes
    steps = Int[]
    append!(steps, fill(stepsizes[1], 6))
    append!(steps, fill(stepsizes[2], 5))
    append!(steps, fill(stepsizes[3], 21))
    append!(steps, fill(stepsizes[4], 48))
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
        PV Generation at last timestep = 0 [kW]
        Actual Power consumption at last timestep = 0 [kW]
        Battery Initial SOC = 0.5
    =#
    global current_j_states = [Initial_PV_Gen, Initial_P_0, Initial_B_SOC*BatterySize, Initial_Curtailment]
    
    Costs = zeros(NumRun);
    AccumulatedCosts = zeros(NumRun);
    Curtailments = zeros(NumRun);
    fmu_temp = zeros(NumRun);
    fmu_indoor_temp = zeros(NumRun);
    Solve_failure = zeros(NumRun);
    AllStates = zeros(NumRun, 4);

    Heat_pump_modes_j = zeros(NumRun);
    HP_H_Usage_j = zeros(NumRun)
    
    fmu_HP_Mode = zeros(NumRun)
    fmu_FCU_Mode = zeros(NumRun)
    fmu_HP_power = zeros(NumRun)
    fmu_Useful_Thermal= zeros(NumRun)
    fmu_HT_PCM_C = zeros(NumRun)
    fmu_HT_PCM_H = zeros(NumRun)
    fmu_system_Mode = zeros(NumRun)
    fmu_FCU = zeros(NumRun)
    fmu_fan = zeros(NumRun)
    fmu_Pump1 = zeros(NumRun)
    fmu_Pump2 = zeros(NumRun)
    fmu_PCM_H_Temp = zeros(NumRun)
    fmu_PCM_C_Temp = zeros(NumRun)
    global TotalCost = 0;

    Actions = zeros(NumRun);
end

############ Data Preparations ############

df = convert_weather(weather, "epw")
Q_total, Q_cd, Q_cv, Q_rh, Q_rc, TC = passive_model(calibration_filepath, df, T_indoor_constant)
pv_cf = Generate_PV(df)
e_load = generate_schedules(f_run, TimeEnd, "complex")


df = combine_dfs(df, Q_total, pv_cf, e_load)


############ Program Execution ############
runtime = @elapsed begin
    starttime = weather.DateTime[TimeStart];
    
    initial_schedule = [0, 0]; # Do nothing while the model is solving the first iteration
    initial_timepoints = [0, 60*stepsize]
    res = FMU.fmu(1, initialstates, initial_schedule, initial_timepoints, model)
    initial_m_states = pandas_to_julia(res)
    
    Actions[1] = 0;   
    AllStates[1, 1] = current_j_states[3]
    AllStates[1, 2] = initial_m_states[end, "PCM_Hot_SOC"]
    AllStates[1, 3] = initial_m_states[end, "PCM_Cold_SOC"]
    
    Curtailments[1] = current_j_states[4]
    fmu_temp[1] = initial_m_states[1, "Ambient Temp"]
    fmu_indoor_temp[1] = initial_m_states[end, "Indoor_Temp"]
    fmu_HP_Mode[1] = initial_m_states[1, "HP Mode"]
    fmu_FCU_Mode[1] = initial_m_states[1, "FCU Mode"]
    fmu_system_Mode[1] = initial_m_states[1, "Exact System Mode"]
    fmu_PCM_H_Temp[1] = initial_m_states[1, "PCM_Hot_Temp"]
    fmu_PCM_C_Temp[1] = initial_m_states[1, "PCM_Cold_Temp"]
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
        timepoints = [Dates.value(weather.DateTime[i]-starttime)/1000, Dates.value(weather.DateTime[i+1]-starttime)/1000] # [s]
        # println(timepoints)

        # input_data = resample_mpc(df, i, steps, stepsize)
        input_data = aggregate_energy(df, i, steps, stepsize)
        println("The input data from the passive model is:")
        println(input_data)
        # Call the MPC function
        currentcost, new_j_states, input_schedule, didnt_solve = OptimizeV4_1(i-1, TC, steps, stepsize, input_data, current_m_states, current_j_states);

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
        fmu_FCU_Mode[i] = new_m_states[1, "FCU Mode"]
        fmu_system_Mode[i] = new_m_states[1, "Exact System Mode"]
        
        fmu_fan[i] = p2e(new_m_states)[1]
        fmu_FCU[i] = p2e(new_m_states)[2]
        fmu_HP_power[i]= p2e(new_m_states)[3]
        fmu_HT_PCM_C[i] = p2e(new_m_states)[4]
        fmu_HT_PCM_H[i] = p2e(new_m_states)[5]
        fmu_Pump1[i] = p2e(new_m_states)[6]
        fmu_Pump2[i] = p2e(new_m_states)[7]
        fmu_Useful_Thermal[i]= p2e(new_m_states)[8]
        fmu_PCM_H_Temp[1] = new_m_states[1, "PCM_Hot_Temp"]
        fmu_PCM_C_Temp[1] = new_m_states[1, "PCM_Cold_Temp"]
        # Is the current iteration of optimization infeasible?
        Solve_failure[i] = didnt_solve
        
        # Store the battery SOC, PCM_H SOC, PCM_C SOC, and 
        
        AllStates[i, 1] = new_j_states[3]
        AllStates[i, 2] = new_m_states[1, "PCM_Hot_SOC"]
        AllStates[i, 3] = new_m_states[1, "PCM_Cold_SOC"]
        Curtailments[i] = new_j_states[4]
        # Store actions
        Actions[i] = input_schedule[1]
        global TotalCost += currentcost;
    end
end

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

# Runtime per iteration (must be less than the MPC run frequency)
rt = runtime/(NumRun-TimeStart)
println("")
println("The runtime is $rt seconds")

println("The total loss of load is $TotalCost kWh")

Total_Inf = sum(Solve_failure)
Inf_pct = 100*Total_Inf/(NumRun-1)

println("There are a total of $Total_Inf infeasible iterations out of $(NumRun-1) iterations, which is $Inf_pct % of the total iterations.")


# EXPORT INTO A CSV FILE (Time series)

############ Save Results to Excel ############
begin
    # Create Excel files
    begin
        # Record the exact time when the result is produced
        timestamp = Dates.format(now(), "yyyy-mm-dd_HH-MM-SS")
        # Create file names and sheet names
        xlsx_file = "MPC$version)_$today_date.$timestamp.xlsx"
        sheetnames = ["TimeSeries Data"]
    end
    
    # Store Result Data
    begin
        column_names = ["DateTime", "Ti (C)", "Ta (C)", "Julia Command", "FMU System Mode", "FMU FCU Mode", "FMU HP Mode", "Curtailment (kWh)",
        "Battery SOC", "Battery Remaining Energy (kWh)", "PCM H SOC", "PCM C SOC", "PCM H Remaining Energy (kWh)", "PCM C Remaining Energy (kWh)",
        "PCM H Temp (C)", "PCM H Remaining Energy (kWh)", ]

        # Collect lists into a tuple or an array
        data_lists = (df[1:NumRun, "DateTime"], fmu_indoor_temp.-273.15, fmu_temp.-273.15, Actions, fmu_system_Mode, fmu_FCU_Mode, fmu_HP_Mode,
        Curtailments, AllStates[:, 1]/BatterySize, AllStates[:, 1], AllStates[:, 2], AllStates[:, 3], AllStates[:, 2].*(PCM_H_Size), AllStates[:, 3].*(PCM_C_Size),
        )  # Fill in with all your lists up to list10
        
        
        Result_Data = DataFrame(
            

            "PCM Hot"= AllStates[:, 2],

            _FMU_Fan_ = fmu_fan,
            _FMU_FCU_ = fmu_FCU,
            _FMU_HP_ = fmu_HP_power,
            _FMU_PCM_H_HeatTransfer_ = fmu_HT_PCM_H,
            _FMU_PCM_C_HeatTransfer_ = fmu_HT_PCM_C,
            _FMU_Pump1_ = fmu_Pump1,
            _FMU_Pump2_ = fmu_Pump2,
            _FMU_Useful_Thermal_ = fmu_Useful_Thermal,
            _Julia_Solver_Failure_ = Solve_failure,
        )
    end
    # Define the path for the Excel file
    excel_file_path = joinpath(folder_path, xlsx_file)

    # Write the DataFrames to different sheets in the same Excel file
    XLSX.openxlsx(excel_file_path, mode="w") do xf

        sheet1 = XLSX.addsheet!(xf, sheetnames[1])
        XLSX.writetable!(sheet1, DataFrames.eachcol(Result_Data), DataFrames.names(Result_Data))
    end
end