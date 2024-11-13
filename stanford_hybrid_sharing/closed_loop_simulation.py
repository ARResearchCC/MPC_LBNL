#!/usr/bin/env python
# coding: utf-8

# In[1]:


from pyfmi import load_fmu
import yaml
import numpy as np
import datetime
import scipy
import matplotlib.pyplot as plt
import pandas as pd
import time
import os
import subprocess

# In[2]:


fmuName='Stanford_Hybrid_System.fmu'


# In[3]:


#Load the FMU
model = load_fmu(fmuName)
fmuNameNoSuffix=fmuName.replace(".fmu","")


# In[4]:


# all time is in seconds
start_time=0
final_time=86400
control_timestep = 21600
data_interval=60


# In[5]:


current_time=start_time
iteration=1

while current_time+control_timestep <= final_time:

    
    
    current_control_timestep_start_time=current_time
    current_control_timestep_final_time=current_time+control_timestep

    
    # this if-else code block is used to save the previous model state for each control timestep
    if iteration==1:
        opts = model.simulate_options()
        opts["ncp"] = (current_control_timestep_final_time-current_control_timestep_start_time)/data_interval
    else:
        state = model.get_fmu_state()
        opts['initialize'] = False
        opts["ncp"] = (current_control_timestep_final_time-current_control_timestep_start_time)/data_interval
        model.set_fmu_state(state)
    
    
    # The "inputs_dict_at_*" variable and the "inputs" variable are used to define the inputs to the FMU model at each control timestep
    inputs_dict_at_start_time={
        'schedule_input':2,
        'mpc_enable':0,
    }
    inputs_dict_at_intermediate_time_1={
        'schedule_input':0,
        'mpc_enable':0,
    }
    inputs_dict_at_intermediate_time_2={
        'schedule_input':3,
        'mpc_enable':0,
    }
    inputs_dict_at_end_time={
        'schedule_input':2,
        'mpc_enable':0,
    }
    inputs = (
        list(inputs_dict_at_start_time), 
        np.array(
            [[current_control_timestep_start_time]+list(inputs_dict_at_start_time.values()),
             [current_control_timestep_start_time+control_timestep*0.2]+list(inputs_dict_at_start_time.values()),
             [current_control_timestep_start_time+control_timestep*0.2]+list(inputs_dict_at_intermediate_time_1.values()),
             [current_control_timestep_start_time+control_timestep*0.6]+list(inputs_dict_at_intermediate_time_1.values()),
             [current_control_timestep_start_time+control_timestep*0.6]+list(inputs_dict_at_intermediate_time_2.values()),
             [current_control_timestep_final_time]+list(inputs_dict_at_intermediate_time_2.values()),
             [current_control_timestep_final_time]+list(inputs_dict_at_end_time.values()),]
        )
    )
    
    # This runs the simulation
    res = model.simulate(start_time=current_control_timestep_start_time, final_time=current_control_timestep_final_time,input=inputs, options=opts)
    
    # This converts the resulting MAT files to csv files
    command =f"python MAT_to_CSV_conversion.py {fmuNameNoSuffix+'_result.mat'} {fmuNameNoSuffix+'_data_points.csv'}"
    subprocess.run(command, shell=True, capture_output=True, text=True)
    
    
    df_current = pd.read_csv(fmuNameNoSuffix+'_result.csv')
    df_current=df_current.drop(['Unnamed: 0'], axis=1)
    
    # The 'output_variables_dict' variable is used to store the simulation outputs at the end of each control timestep
    # If you want to apply external control to this simulation, you could find a way to use 
    # the 'output_variables_dict' variable to control the 'inputs' variable
    output_variables_dict=df_current.iloc[-1].to_dict()
    
    if iteration==1:
        df = pd.DataFrame()

    df=pd.concat([df[:-1], df_current], ignore_index=True)    
    
    iteration+=1
    current_time += control_timestep


# In[6]:


try:
    os.remove(fmuNameNoSuffix+'_result.mat')
except:
    pass
try:
    os.remove(fmuNameNoSuffix+'_result.csv')
except:
    pass

try:
    os.remove(fmuNameNoSuffix+'_log.txt')
except:
    pass


# In[7]:


df.to_csv(fmuNameNoSuffix+'_result.csv',index=False)


# In[ ]:





# # Plotting

# In[8]:


fig,ax=plt.subplots(nrows=1,ncols=1,figsize=(10,7))

ax.plot(df['time'],df['schedule_input'], label="schedule_input",color="red")



ax.legend(fontsize=15)
ax.set_title(f'schedule_input Plot',fontsize=20)
ax.set_xlabel('Time',fontsize=18)
ax.set_ylabel('schedule_input',fontsize=18)


# In[ ]:




