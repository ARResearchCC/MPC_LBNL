# from pyfmi import load_fmu
import numpy as np
import pandas as pd

'''
This function is being called by Julia MPC (MPC_S V4.jl) through PyCall. The fmu function 
will take the following inputs:
1. run iteration (currently unused)
2. modelica states (hot_PCM.T_start, cold_PCM.T_start, RoomModel.roo.T_start etc.) initialized
    ONLY at iteration = 1. At other iterations, the fmu will get the end states of last iteration.
3. commands (0, 1, 2, 3) optimized by Julia MPC to operate the system in Modelica through FMU
4. time points which include a time series array from the starting time point to the final
    time point of the current run iteration.
5. the fmu model
This function will call the FMU which has the Modelica model, and given the current states of
the Modelica model (passed as input of the FMU function at i=1 or remembered itself in the model from last iteration 
after i=1) and the operation commands, execute the Modelica simulation which will return the "ground truth" states 
including electricity consumption from the HP-PCM system, accurate SOC of the PCMs, and indoor temperature etc. 
These states are called M_States, which will be fed into the Julia MPC in the next iteration.

'''

def fmu(iteration, states, commands, time_points, model):

    # Define the inputs for the simulation
    input_names = ['schedule_input', 'mpc_enable']
    mpc_enable_values = np.ones(len(time_points))  # Example values for 'mpc_enable'
    
    # Combine time points with input values
    input_values = np.column_stack((time_points, commands, mpc_enable_values))
    
    # Reformat inputs for simulation
    input_data = (input_names, input_values)

    print(input_data)
    opts = model.simulate_options()
    # Initialize Modelica model states at iteration = 1 or use last saved states from last iteration
    try:
        if iteration == 1:

            # Set up the experiment at iteration = 1
            model.set('hot_PCM.T_start', states[0])
            model.set('cold_PCM.T_start', states[1])
            model.set('RoomModel.roo.T_start', states[2])
            model.set('hot_PCM.Tes_nominal', states[3]*3600000)
            model.set('cold_PCM.Tes_nominal', states[4]*3600000)
            model.set('fanCoilUnit_hysteresis_air_control.zone_temp_cooling_setpoint', states[5])
            model.set('fanCoilUnit_hysteresis_air_control.zone_temp_heating_setpoint', states[6])
            model.set('fanCoilUnit_hysteresis_air_control.zone_temp_setpoint_delta', states[7])
            model.set("hybrid_2024_controller_new.modular_SOC_Manager1.hysteresis[1].uLow", 0)
            model.set("hybrid_2024_controller_new.modular_SOC_Manager1.hysteresis[1].uHigh", 0.01)
            model.set("hybrid_2024_controller_new.modular_SOC_Manager1.hysteresis[2].uLow", 0)
            model.set("hybrid_2024_controller_new.modular_SOC_Manager1.hysteresis[2].uHigh", 0.01)

            opts['initialize'] = True
            opts["ncp"] = len(time_points)
            
        else:
            
            # Retrieve fmu model ending state from last iteration if iteration > 1
            state = model.get_fmu_state()
            opts['initialize'] = False
            opts["ncp"] = len(time_points)
            model.set_fmu_state(state)

        # Run the simulation from start_time to final_time with specified inputs and options
        res = model.simulate(start_time=time_points[0], final_time=time_points[-1], input=input_data, options=opts)
    except Exception as e:
        print(f"An error occurred: {e}")
        raise

    # Collect output
    time = res['time']
    data = {
        'Time': time, 
        'Fan Coil Fan Electric Power': res['fan_coil_fan_electric_power'], 
        'HP Electric Power': res['heat_pump_electric_power'],
        'HP Mode': res['heat_pump_mode'], 
        'HP ON OFF': res['heat_pump_on_off'], 
        'HP Useful Thermal Power': res['heat_pump_useful_thermal_power'],
        'Ambient Temp': res['outdoor_air_temperature'], 
        'Pump 1 Electric Power': res['pump_1_electric_power'], 
        'Pump 2 Electric Power': res['pump_2_electric_power'], 
        'Indoor_Temp': res['room_temperature'], 
        'PCM_Cold_SOC': res['soc_cold_pcm'], 
        'PCM_Hot_SOC': res['soc_hot_pcm'], 
        'PCM_Cold_Temp': res['temp_cold_pcm'], 
        'PCM_Hot_Temp': res['temp_hot_pcm'],
        'Heat Transfer Cold PCM': res['heat_transfer_cold_pcm'], 
        'Heat Transfer Hot PCM': res['heat_transfer_hot_pcm'], 
        'FCU Heat Delivered': res['fcu_heat_delivered'], 
        'Exact System Mode': res['exact_system_mode'], 
        'FCU Entering Water Temperature ': res['fcu_entering_water_temperature'], 
        'FCU Leaving Water Temperature ': res['fcu_leaving_water_temperature'],
        'Fan Coil Fan Mass Flow': res['fan_coil_fan_mass_flow'], 
        'Heat Pump Supply Water Temperature': res['heat_pump_supply_water_temperature'],
        'Heat Pump Return Water Temperature': res['heat_pump_return_water_temperature'], 
        'Heat Pump Compressor Speed': res['heat_pump_compressor_speed'],
        'FCU Mode': res['FMUMode']
    }

    # This is the M_States output used for Julia at next iteration:
    df = pd.DataFrame(data)
    return df


'''
current_m_states = [273.15+48, 273.15+11, 273.15+22, 11, 7, 273.15+22, 273.15+20, 0.5, 0.5, 0.5]

commands = [[1, 1], [0, 0]]
times = [[0, 3600], [3600, 7200]]
for i in [1, 2]:
    df = fmu(i, current_m_states, commands[i-1], times[i-1])
    print(df)
    print(df['PCM_Hot_Temp'])
'''