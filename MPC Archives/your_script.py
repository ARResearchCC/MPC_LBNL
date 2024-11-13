from pyfmi import load_fmu
import numpy as np
import matplotlib.pyplot as plt
import csv
import pandas as pd

def fmu(a):
    # Define the path to the FMU file
    fmuName = "C:/Users/Fred/Desktop/PyFMI/Stanford_Hybrid_System.fmu"
    
    # Load the FMU model
    model = load_fmu(fmuName)
    
    # Define the inputs for the simulation
    input_names = ['schedule_input', 'mpc_enable']
    
    # Define time points and corresponding input values
    time_points = [0.0, 300.0, 600.0]  # Example time points
    schedule_input_values = [3, 3, 3]  # Example values for 'schedule_input'
    # schedule_input_values = [0, 0, 0]
    mpc_enable_values = [1, 1, 1]  # Example values for 'mpc_enable'
    
    # Combine time points with input values
    input_values = np.column_stack((time_points, schedule_input_values, mpc_enable_values))
    
    # Reformat inputs for simulation
    input_data = (input_names, input_values)
    
    # Set up the experiment
    # model.setup_experiment(start_time=0, stop_time=600)
    
    model.set('hot_PCM.T_start', 273.15 + 55)
    model.set('cold_PCM.T_start', 273.15 + 10)
    model.set('RoomModel.roo.T_start', 273.15 + 21)
    # Retrieve simulation options and current FMU state
    opts = model.simulate_options()
    
    # Modify simulation options
    opts['initialize'] = True
    opts["ncp"] = 100  # Number of communication points
    
    # Run the simulation from start_time to final_time with specified inputs and options
    res = model.simulate(start_time=0, final_time=600, input=input_data, options=opts)
    
    time = res['time']

    data = {'Time': time, 'Fan Coil Fan Electric Power': res['fan_coil_fan_electric_power'], 'Heat Pump Electric Power': res['heat_pump_electric_power'],
            'Heat Pump Mode': res['heat_pump_mode'], 'Heat Pump ON OFF': res['heat_pump_on_off'], 'Heat Pump Useful Thermal Power': res['heat_pump_useful_thermal_power'],
            'MPC Enable': res['mpc_enable'], 'Ambient Temp': res['outdoor_air_temperature'], 'Pump 1 Electric Power': res['pump_1_electric_power'], 
            'Pump 2 Electric Power': res['pump_2_electric_power'], 'Indoor Temp': res['room_temperature'], 'Schedule Input': res['schedule_input'],
            'PCM Cold SOC': res['soc_cold_pcm'], 'PCM Hot SOC': res['soc_hot_pcm']}
    df = pd.DataFrame(data)

    return df

# Call the function with an example argument
df = fmu(1)

print(df)


'''
# Write the keys to a CSV file
csv_file_path = "C:/Users/Fred/Desktop/PyFMI/fmu_keys.csv"
with open(csv_file_path, mode='w', newline='') as file:
    writer = csv.writer(file)
    for key in keys_list:
        writer.writerow([key])

# Print the keys list
print("Keys in res:", keys_list)
print(f"Keys have been written to {csv_file_path}")
'''
