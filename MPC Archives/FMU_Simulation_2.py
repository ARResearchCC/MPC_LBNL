from pyfmi import load_fmu
import numpy as np
import pandas as pd

'''
This function is being called by Julia MPC (MPC_S V2.jl) through PyCall. The fmu function 
will take the following inputs:
1. run iteration (currently unused)
2. modelica states (hot_PCM.T_start, cold_PCM.T_start, RoomModel.roo.T_start) either initialized
    at iteration = 1 or the output from the last iteration of the FMU function.
3. commands (0, 1, 2, 3) optimized by Julia MPC to operate the system in Modelica through FMU
4. time points which include a time series array from the starting time point to the final
    time point of the current run iteration.
This function will call the FMU which has the Modelica model, and given the current states of
the Modelica model (passed as input of the FMU function) and the operation commands, execute the 
Modelica simulation which will return the "ground truth" states including electricity consumption 
from the HP-PCM system, accurate SOC of the thermal batteries, and indoor temperature etc. 
The additional output parameters from the FMU compared to V1 are PCM Hot Temperature and PCM Cold Temperature.
These states are called M_States, which will be fed into the Julia MPC in the next iteration, as well
as the FMU in the next iteration.

'''

def fmu(iteration, states, commands, time_points):

    # Define the path to the FMU file
    # fmuName = "C:/Users/Fred/Desktop/PyFMI/Stanford_Hybrid_System.fmu"
    fmuName = "C:/Users/Fred/Desktop/PyFMI/Stanford_Hybrid_System_2.fmu"

    # Load the FMU model
    model = load_fmu(fmuName)
    
    # Define the inputs for the simulation
    input_names = ['schedule_input', 'mpc_enable']
    
    # Define time points and corresponding input values
    # time_points = [0.0, 30.0, 60.0]  # Example time points
    # schedule_input_values = [1, 0, 1]  # Example values for 'schedule_input'
    mpc_enable_values = np.ones(len(time_points))  # Example values for 'mpc_enable'
    
    # Combine time points with input values
    input_values = np.column_stack((time_points, commands, mpc_enable_values))
    
    # Reformat inputs for simulation
    input_data = (input_names, input_values)

    # Set up the experiment
    model.setup_experiment(start_time=time_points[0], stop_time=time_points[-1])
    
    # Initialize Modelica model states at iteration = 1 or use last saved states from last iteration
    try:
        opts = model.simulate_options()
        opts["ncp"] = len(time_points)
        opts['initialize'] = True
        model.set('hot_PCM.T_start', states[0])
        model.set('cold_PCM.T_start', states[1])
        model.set('RoomModel.roo.T_start', states[2])

        print("About to start FMU simulation")
        # print("Start time:", time_points[0])
        # print("Final time:", time_points[-1])
        # print("Input data is: ")
        # print(input_data)
        # print("Simulation options:", opts)
        
        # Run the simulation from start_time to final_time with specified inputs and options
        res = model.simulate(start_time=time_points[0], final_time=time_points[-1], input=input_data, options=opts)
    except Exception as e:
        print(f"An error occurred: {e}")
        raise
    
    # Collect output
    time = res['time']

    data = {'Time': time, 'Fan Coil Fan Electric Power': res['fan_coil_fan_electric_power'], 'HP Electric Power': res['heat_pump_electric_power'],
            'HP Mode': res['heat_pump_mode'], 'HP ON OFF': res['heat_pump_on_off'], 'HP Useful Thermal Power': res['heat_pump_useful_thermal_power'],
            'MPC Enable': res['mpc_enable'], 'Ambient Temp': res['outdoor_air_temperature'], 'Pump 1 Electric Power': res['pump_1_electric_power'], 
            'Pump 2 Electric Power': res['pump_2_electric_power'], 'Indoor_Temp': res['room_temperature'], 'Schedule Input': res['schedule_input'],
            'PCM_Cold_SOC': res['soc_cold_pcm'], 'PCM_Hot_SOC': res['soc_hot_pcm'], 'PCM_Cold_Temp': res['temp_cold_pcm'], 'PCM_Hot_Temp': res['temp_hot_pcm']}
    df = pd.DataFrame(data)

    return df
