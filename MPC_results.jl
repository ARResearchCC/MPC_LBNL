# Model Predictive Control Result
# Complete Date: 11/18/2024
# ACE Microgrid Project --- CEE Atmosphere and Energy --- Stanford University

# This code is an intellectual property of Yuanbei "Fred" Fan.
# Dear user, this is Fred. If there are any questions with regards to the code, please do not feel hesistant to reach out to me at yf1098@stanford.edu.

############ Initialize Tools ############
begin
    using Pkg

    # Please unmute the next two lines if you are running this code for the first time and don't have all the packages installed.
    # Pkg.add("JuMP", "CSV", "DataFrames", "PlotlyJS", "Dates", "XLSX", "FileIO", "Base", "Random", "Statistics")
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
end

############ Load Data ############
using XLSX

# Load the Excel file
file_path ="C:\\Users\\Fred\\Desktop\\PyFMI\\MPC_S_V5.0\\MPC(5.0)(2024-11-21_10-08-11).xlsx"

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

Continuous_Data = dfs["Continuous Data"]
Discrete_Data = dfs["Discrete Data"]

# Initialize Ideal_Mode
Ideal_Mode = zeros(size(Discrete_Data, 1))

# Loop through each row of Discrete_Data
for i in 1:size(Discrete_Data, 1)
    julia_command = Discrete_Data[i, "MPC Command"]
    fmu_mode = Discrete_Data[i, "FMU FCU Mode"]

    if julia_command == 0 && fmu_mode == -1
        Ideal_Mode[i] = 10
    elseif julia_command == 0 && fmu_mode == 0
        Ideal_Mode[i] = 0
    elseif julia_command == 0 && fmu_mode == 1
        Ideal_Mode[i] = 9
    elseif julia_command == 1 && fmu_mode == -1
        Ideal_Mode[i] = 10
    elseif julia_command == 1 && fmu_mode == 0
        Ideal_Mode[i] = 1
    elseif julia_command == 1 && fmu_mode == 1
        Ideal_Mode[i] = 3
    elseif julia_command == 2 && fmu_mode == -1
        Ideal_Mode[i] = 7
    elseif julia_command == 2 && fmu_mode == 0
        Ideal_Mode[i] = 5
    elseif julia_command == 2 && fmu_mode == 1
        Ideal_Mode[i] = 9
    elseif julia_command == 3 && fmu_mode == -1
        Ideal_Mode[i] = 6
    elseif julia_command == 3 && fmu_mode == 0
        Ideal_Mode[i] = 0
    elseif julia_command == 3 && fmu_mode == 1
        Ideal_Mode[i] = 2
    end
end


# Plot Accumulated loss of load
trace1 = scatter(
    x = Discrete_Data[:, "DateTime"],  
    y = Ideal_Mode,
    name="FMU System Mode in theory"
)

trace2 = scatter(
    x = Discrete_Data[:, "DateTime"],  
    y = Discrete_Data[:, "FMU System Mode"],
    name="FMU System Mode in reality"
)

trace3 = scatter(
    x = Discrete_Data[:, "DateTime"],  
    y = Discrete_Data[:, "FMU System Mode"] .- Ideal_Mode,
    name="FMU System Mode in reality"
)




trace4 = scatter(
    x = Discrete_Data[:, "DateTime"],  
    y = Discrete_Data[:, "PCM H Temp (C)"],
    name="PCM H Temp (C)"
)

trace5 = scatter(
    x = Discrete_Data[:, "DateTime"],  
    y = Discrete_Data[:, "PCM C Temp (C)"],
    name="PCM C Temp (C)"
)

trace6 = scatter(
    x = Discrete_Data[:, "DateTime"],  
    y = Discrete_Data[:, "PCM H SOC"],
    name="PCM H SOC"
)

trace7 = scatter(
    x = Discrete_Data[:, "DateTime"],  
    y = Discrete_Data[:, "MPC Command"],
    name="Julia Command"
)

trace8 = scatter(
    x = Discrete_Data[:, "DateTime"],  
    y = Discrete_Data[:, "FCU Ave Heat Power (kW)"],
    name="FCU Ave Heat Power (kW)"
)
trace9 = scatter(
    x = Discrete_Data[:, "DateTime"],  
    y = Discrete_Data[:, "Ti (C)"],
    name="indoor temp"
)


p_example = plot([trace8, trace9, trace2, trace7], Layout(title="FCU Ave Heat Power (kW)", xaxis_title="Time", yaxis_title="FCU Ave Heat Power (kW)"))
display(p_example)

p_modes = plot([trace1, trace2], Layout(title="Theoretical and Real FMU System Mode", xaxis_title="Time", yaxis_title="Mode"))
display(p_modes)

p_modes_diff = plot([trace3, trace6], Layout(title="Mode Difference", xaxis_title="Time", yaxis_title="Mode"))
display(p_modes_diff)