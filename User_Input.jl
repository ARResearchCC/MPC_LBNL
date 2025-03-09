# User Input File

Python_Environment = "C:\\Users\\Fred\\anaconda3\\python.exe"
Gurobi_License = "C:\\Users\\Fred\\gurobi.lic"

fmu_Name = "C:\\Users\\Fred\\Desktop\\PyFMI\\Stanford_Hybrid_System.fmu"
FMU_Simulation_File = "FMU_Simulation_4"

#=
version = 5.0
code_name = "MPC_S_V$version"
=#

version = 1.0
code_name = "Rule_Based_V$version"

# Load the weather file. This weather file is a historical weather file with 5 minutes step size. 
# This weather file is treated as the ground truth, but simulated forecast errors are added later.
weather_file_name = "C:\\Users\\Fred\\Desktop\\Microgrid_Project\\Data Cleaning\\updated_epw_weather_nonleap_5mins.csv"

calibration_file_path = "C:\\Users\\Fred\\Desktop\\PyFMI\\Calibration_Model_Input.xlsx"

Main_Folder_Path = "C:\\Users\\Fred\\Desktop\\PyFMI"

TimeStart = 1;
f_run = 12 # run frequency [times/hour]
NumRun = 12*18 # 1 week
δt = 1/f_run # [hr]
# NumRun = (TimeEnd-TimeStart+1) - Opt_Horizon + 1; # This is the max number of runs allowed for a full year of data

# Define the sampling intervals
# stepsizes = [5, 30, 60, 120] # minutes
# stepnums = [6, 5, 21, 48]
# stepsizes = [60] # minutes
# stepnums = [72]
