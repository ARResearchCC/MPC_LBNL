# MPC Algorithm
#=
This function is being called by Julia MPC (MPC_S V4.jl) as the actual MPC algorithm. The MPC
will take the following inputs:
1. run iteration
2. stepsize (5 mins)
3. weather file (1 week long starting at the current timestep) which includes
    parameters such as ambient temperature, wind speed, solar irradiance etc.
4. the electrical load schedules (plug loads and lighting)
5. M_states (either initiated at iteration = 1, or taken from the previous iteration of the Modelica
model) including the electricity consumption from the HP-PCM system, accurate SOC of the thermal 
batteries, and indoor temperature etc.
6. J_States (either initiated at iteration = 1, or taken from the previous iteration of the MPC
model) including the actual solar generation, lighting & plug load consumption from last iteration,
and the accurate starting battery SOC from the previous iteration.

The J_States basically include the parameters that are not updated from Modelica.

This function first compute the accurate starting battery SOC at current iteration based on the 
M_States and J_States, and use the M_States and J_States (both assumed as accurate ground truth)
parameters to initialize the current MPC model. The MPC model will use a reduced order model as 
basis to make optimal decisions (operational commands for HP-PCM system) which will be used as input
for the current iteration of Modelica simulation in the FMU. It will also return the J_States which
will be used as input for the next iteration of the MPC.
=#
using Gurobi
function OptimizeV4(iteration, δt, Weather, Schedules, M_States, J_States)
    ########## Instructions  ##########
    # current optimization model output for modelica (0 to 3)
    ENV["GRB_LICENSE_FILE"] = "C:\\Users\\Fred\\gurobi.lic"
    # ENV["GRB_LICENSE_FILE"] = "C:\\Users\\Fred\\.julia\\environments\\v1.10\\gurobi.lic" 
    ########## Data Preparations  ##########  
    begin
        # Set timesteps 
        NumTime = size(Weather)[1] # 168 * 12 intervals for now

        # Add "fake" forecast errors to the ambient temperature forecast
        Weather.Temperature = Add_Uncertainty_Temp(Weather.Temperature);

        Weather.CFM = zeros(NumTime)
        Weather.RadHeatGain = zeros(NumTime)
        Weather.RadCool = zeros(NumTime)
        Weather.PV = zeros(NumTime)

        for i = 1:NumTime
            Weather.CFM[i] = Al*sqrt(Cs*(T_indoor_constant + 273.15 - T_d) + Cw*wind_speed[i]^2) * 2.11888 # [ft^3/min] 
            Weather.RadHeatGain[i] = calculate_solarheatgain(Weather[i,:datetime], Weather[i,:DNI], Weather[i,:DHI], Weather[i,:GHI]) # [BTU/hr]
            day_status = IsDay(Weather.SolarTime[i], Weather.datetime[i]) 
            Weather.RadCool[i] = Q_radiativeCooling(T_indoor_constant, Weather.Temperature[i], Weather[i,:"Dew Point"], Weather.SolarTime[i], day_status)
            Weather.PV[i] = calculate_PV(Weather[i,:datetime], Weather[i,:DNI], Weather[i,:DHI], Weather[i,:GHI], 180-Berg_tilt, module_surface_tilt, Weather[i,:celltemp])
        end

        # Declare time series variables
        TemperatureAmbient = Weather.Temperature;
        PercentLighting = Schedules[:,1];
        PercentPlug = Schedules[:,2];
        PercentOccupied = Schedules[:,3];
        PVGeneration = Weather.PV;
        RadHeatGain = Weather.RadHeatGain;
        RadCooling = Weather.RadCool;
        CFM = Weather.CFM;

        # Length of the M_States dataframe
        msize = size(M_States)[1];

        # Initialize the starting values with information from J_States and M_States
        InStoragePCM_C_1 = min(M_States[msize, :"PCM_Cold_SOC"], 1.05) # [1]
        InStoragePCM_H_1 = min(M_States[msize, :"PCM_Hot_SOC"], 1.05) # [1]
        InStorageBattery_1 = J_States[3] # [kWh]
        TemperatureIndoor_1 = (M_States[msize, :"Indoor_Temp"] - 273.15) * (9/5) + 32 # [°F]
        
        # These two variables mean that weather the system is in heating mode or cooling mode at the start of the current iteration.
        heating_0 = 0
        cooling_0 = 0

        if iteration > 1
            println("M_States")
            print(M_States)
            println()
            time_intervals = zeros(msize-1); # [hr]
            
            heating_0 = 0
            cooling_0 = 0
            
            if M_States[msize, :"HP ON OFF"] == 1 && M_States[msize, :"HP Mode"] == 1 
                heating_0 = 1
            end
            if M_States[msize, :"HP ON OFF"] == 1 && M_States[msize, :"HP Mode"] == 0
                cooling_0 = 1
            end

            for i = 1:msize-1
                time_intervals[i] = (M_States[i+1, :"Time"] - M_States[i, :"Time"])/3600 # [hr]
            end
            println("time intervals")
            print(time_intervals)
            println()
            println("E_Pump_1_0")
            E_Pump_1_0 = sum(((M_States[i+1, :"Pump 1 Electric Power"] + M_States[i, :"Pump 1 Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            print(E_Pump_1_0)
            println()
            E_Pump_2_0 = sum(((M_States[i+1, :"Pump 2 Electric Power"] + M_States[i, :"Pump 2 Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            E_HP_0 = sum(((M_States[i+1, :"HP Electric Power"] + M_States[i, :"HP Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            Thermal_Power_Delivered_0 = sum(((M_States[i+1, :"HP Useful Thermal Power"] + M_States[i, :"HP Useful Thermal Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            E_Fan_0 = sum(((M_States[i+1, :"Fan Coil Fan Electric Power"] + M_States[i, :"Fan Coil Fan Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            
            # Julia States
            PVGen_0 = J_States[1] # PV Generation at last timestep (Actual)
            E_other_0 = J_States[2] # [kW] Other Electrical Load (from schedule) at last timestep (Actual)
            B_SOC_0 = J_States[3] # Battery SOC at last timestep (Actual)
            
            E_Modelica_0 = E_Pump_1_0 + E_Pump_2_0 + E_HP_0 + E_Fan_0 # [kWh] Actual Electrical Load from HP system at last timestep
            Total_Load_0 = E_Modelica_0 + E_other_0 # [kWh] Actual total Electrical Load at last timestep
            
            # Small Optimization Model to determine initial battery SOC using J_States
            begin
                println("Start solving small model")

                m_0 = Model(Gurobi.Optimizer)
                # m_0 = Model(Clp.Optimizer)
                # define unit
                @variable(m_0, InStorageBattery_past[1:2] >= 0); # [kWh]
                @variable(m_0, Curtailment_past >= 0); # [kW]
                @variable(m_0, PV2B_past >= 0); # [kW]
                @variable(m_0, PV2H_past >= 0); # [kW]
                @variable(m_0, B2H_past >= 0); # [kW]
                @objective(m_0, Min, Curtailment_past); # [kW]

                @constraint(m_0, PVGen_0 == PV2B_past + PV2H_past + Curtailment_past) # [kW]
                @constraint(m_0, InStorageBattery_past[2] == InStorageBattery_past[1] * δt * (1 - BatteryLoss) +  δt * (PV2B_past * η - B2H_past)); # [kWh]
                @constraint(m_0, Total_Load_0 == PV2H_past * η_PVIV * δt + B2H_past * η * δt) # [kWh]

                # Battery discharging constraint, node at battery 
                @constraint(m_0, InStorageBattery_past[1] == B_SOC_0); # [kWh]

                # Battery discharging constraint, node at battery 
                @constraint(m_0, δt * B2H_past <= InStorageBattery_past[1]); # [kWh]
                
                # Battery power inverter constraint, node at battery (inverter power constraint) 
                @constraint(m_0, B2H_past + PV2B_past <= InverterSize); # [kW]
                
                # Battery storage size constraint, node at battery
                @constraint(m_0, [t=1:2], InStorageBattery_past[t] <= BatterySize); # [kWh]
                
                # Battery storage max discharge constraint, node at battery, always at least 20% full 
                @constraint(m_0, [t=1:2], InStorageBattery_past[t] >= BatterySize * (1-MaxDischarge)); # [kWh]

                optimize!(m_0);
                
                # Actual Battery SOC to start with
                InStorageBattery_1 = value.(InStorageBattery_past[2]) # [kWh]
                println("Finished solving small model")
                #=
                println()
                println("Battery SOC is:")
                print(InStorageBattery_1)
                println()
                println("Hot PCM SOC is:")
                print(InStoragePCM_H_1)
                println()
                println("Cold PCM SOC is:")
                print(InStoragePCM_C_1)
                println()
                println("Indoor Temperature is:")
                print(TemperatureIndoor_1)
                println()
                =#
            end
        end    
    end
    println("Start solving big model")  
    ########## Declare model  ##########
    begin
        # Define the model name and solver. In this case, model name is "m"
        # m = Model(Clp.Optimizer)
        # m = Model(Ipopt.Optimizer)
        begin
            # Set path to license (for those using Gurobi)
            m = Model(Gurobi.Optimizer)
        end
    end
    ######## Decision variables ########
    begin
        @variable(m, PV2H[1:NumTime] >= 0); # [kW] electrical power transfer from PV to home (Berg)

        @variable(m, PV2G[1:NumTime] >= 0); # [kW] electrical power transfer from PV to ground (curtailment)

        @variable(m, PV2B[1:NumTime] >= 0); # [kW] electrical power transfer from PV to battery

        @variable(m, B2H[1:NumTime] >= 0); # [kW] electrical power transfer from battery to home (Berg)
        
        @variable(m, HP_state[1:3, 1:NumTime], Bin) # stand alone = 1, heat space = 2, cool space = 3

        @variable(m, heating[1:NumTime], Bin);
        @variable(m, heating_start[1:NumTime], Bin);
        @variable(m, heating_end[1:NumTime], Bin);
        @variable(m, heating_start_binary[1:NumTime], Bin)
        @variable(m, heating_end_binary[1:NumTime], Bin)
        
        @variable(m, cooling[1:NumTime], Bin);
        @variable(m, cooling_start[1:NumTime], Bin);
        @variable(m, cooling_end[1:NumTime], Bin);
        @variable(m, cooling_start_binary[1:NumTime], Bin)
        @variable(m, cooling_end_binary[1:NumTime], Bin)

        @variable(m, PCM_state[1:5, 1:NumTime], Bin) # Stand alone = 1, charge hot = 2, charge cold = 3, discharge hot = 4, discharge cold = 5

        @variable(m, PCM_SOC[1:2, 1:NumTime] >= 0); # Hot = 1, Cold = 2

        @variable(m, PCM_H_charge_rate[1:NumTime] >= 0); # [BTU/hr] thermal power transfer from HP heating to PCM hot storage

        @variable(m, PCM_C_charge_rate[1:NumTime] >= 0); # [BTU/hr] thermal power transfer from HP cooling to PCM cold storage

        @variable(m, InStorageBattery[1:NumTime] >= 0); # [kWh] Battery Remaining Charge

        @variable(m, TemperatureIndoor[1:NumTime] >= 0); # [°F] Indoor Air Temperature

        @variable(m, G2H[1:NumTime] >= 0); # [kW] Loss of Load
    end
    ############ Objective Functions #############
    begin
        # Set single objective for minimizing annual total cost

        # Penalty for loss of load [$]
        @expression(m, penalty_l, M * δt * sum(G2H[t] for t = 1:NumTime))

        # Total Cost over Optimization Horizon [$]
        # @objective(m, Min, short_run_marginal_cost + penalty_c + penalty_l);
        @objective(m, Min, δt * sum(G2H[t] for t = 1:NumTime));
    end
    ############# Expressions ############
    begin
        # DeltaTemp 
        @expression(m, TempDelta[t=1:NumTime], TemperatureAmbient[t] - TemperatureIndoor[t]); # [°F]

        # Electricity usage from lighting 
        @expression(m, E_Lighting[t=1:NumTime], PeakLighting * PercentLighting[t]); # [kW]

        # Electricity usage from plugs
        @expression(m, E_Plugs[t=1:NumTime], PeakPlugLoad * PercentPlug[t]); # [kW]

        # Electricity usage from Pump
        @expression(m, E_Pumps[t=1:NumTime], sum(PCM_state[ω, t] for ω = 2:5) * P_Pumps); # [kW]

        # Total electricity usage
        @expression(m, E_total[t=1:NumTime], E_Lighting[t] + E_Plugs[t] + P_Controls + E_Pumps[t] + P_AHU); # [kW]

        # Calculate Ventilation
        @expression(m, CFMVen[t=1:NumTime], min(Rp * PercentOccupied[t] * MaxOccupancy + Ra * Area, PercentOccupied[t] * MaxOccupancy * Ventilation)); # [ft^3/min]

        # Heat gain from lighting
        @expression(m, Q_Lighting[t=1:NumTime], PeakLighting * PercentLighting[t] * 3412.14); # [BTU/hr]

        # Heat gain from plugs
        @expression(m, Q_Plugs[t=1:NumTime], PeakPlugLoad * PercentPlug[t] * 3412.14); # [BTU/hr]

        # Heat gain from occupancy
        @expression(m, Q_Occupancy[t=1:NumTime], MaxOccupancy * PercentOccupied[t] * TotalPersonHeat); # [BTU/hr]

        # Heat gain from infiltration
        # 1.08 = specific heat capacity of air at STP:0.24 [BTU/(lb*°F)] * air density at STP:0.075 [lb/ft^3] * 60 [min/hr] 
        @expression(m, Q_Infiltration[t=1:NumTime], 1.08 * TempDelta[t] * CFM[t]); # [BTU/hr]

        # Heat gain from ventilation
        @expression(m, Q_Ventilation[t=1:NumTime], 1.08 * TempDelta[t] * CFMVen[t]); # [BTU/hr]

        # Heat gain from lighting, plugs, occupancy, ventilation, infiltration
        @expression(m, Q_Others[t=1:NumTime], Q_Infiltration[t] + Q_Ventilation[t] + Q_Occupancy[t] + Q_Lighting[t] + Q_Plugs[t]); # [BTU/hr]

        # Heat gain through structural evenlope
        @expression(m, Q_Envelope[t=1:NumTime], UA * TempDelta[t]); # [BTU/hr]

        # Heat gain from solar radiation
        @expression(m, Q_Rad[t=1:NumTime], SHGC * RadHeatGain[t]); # [BTU/hr]

        # Radiative cooling
        @expression(m, Q_RadCool[t=1:NumTime], RCC * RadCooling[t]); # [BTU/hr]

        # Electric Power consumed by Heat Pump (either heating or cooling)
        @expression(m, H2HP[t=1:NumTime], HP_state[2, t] * HP_power_H); # [kW]
        @expression(m, H2C[t=1:NumTime], HP_state[3, t] * HP_power_C); # [kW]

        # Charging and discharging of PCM storages
        @expression(m, PCM_H_Charge[t=1:NumTime], PCM_state[2, t] * PCM_H_charge_rate[t]); # [BTU/hr]
        @expression(m, PCM_C_Charge[t=1:NumTime], PCM_state[3, t] * PCM_C_charge_rate[t]); # [BTU/hr]
        @expression(m, PCM_H_Discharge[t=1:NumTime], PCM_state[4, t] * PCM_H_discharge_rate); # [BTU/hr]
        @expression(m, PCM_C_Discharge[t=1:NumTime], PCM_state[5, t] * PCM_C_discharge_rate); # [BTU/hr]

        # Heat or "Coolth" delivered to the room (directly from Heat Pump or PCM discharge)
        @expression(m, HP2R[t=1:NumTime], H2HP[t] * COP_H * 3412.14 - PCM_H_Charge[t]); # [BTU/hr]
        @expression(m, C2R[t=1:NumTime], H2C[t] * COP_C * 3412.14 - PCM_C_Charge[t]); # [BTU/hr]


        # Heat or "Coolth" delivered to the room (directly from Heat Pump or PCM discharge)
        @expression(m, Heat_delivered[t=1:NumTime], HP2R[t] + PCM_H_Discharge[t]); # [BTU/hr]
        @expression(m, Cool_delivered[t=1:NumTime], C2R[t] + PCM_C_Discharge[t]); # [BTU/hr]
     

    end
    ############# Constraints ############
    begin
        # Initialization Constraints with States from last timestep
        begin
            # Set point temperature range constraints
            @constraint(m, TemperatureIndoor[1] == TemperatureIndoor_1); # [°F]

            # Battery storage initialization constraint
            @constraint(m, InStorageBattery[1] == InStorageBattery_1); # [kWh] 
        end

        # Heat pump operational constraint
        @constraint(m, [t=1:NumTime], sum(HP_state[:, t]) == 1); # Heat Pump can only be in one mode at a time

        # Thermal storage operational constraint
        @constraint(m, [t=1:NumTime], sum(PCM_state[:, t]) == 1); # PCM can only be in one mode at a time

        # Internal temperature balance evolution constraint
        @constraint(m, [t=1:NumTime-1], TemperatureIndoor[t+1] == TemperatureIndoor[t] + 
        (δt/TC)*(Q_Others[t] + Q_Envelope[t] + Q_Rad[t] - Q_RadCool[t] + Heat_delivered[t] - Cool_delivered[t])); # [°F]
        
        # PV energy balance constraint, node at PV
        @constraint(m, [t=1:NumTime], PVGeneration[t] * PVSize == PV2B[t] + PV2H[t] + PV2G[t]); # [kW]
    
        # House electricity load constraint, node at house, battery efficiency modeled
        @constraint(m, [t=1:NumTime], E_total[t] + H2HP[t] + H2C[t] == PV2H[t] * η_PVIV + B2H[t] * η + G2H[t]); # [kW]

        # Battery storage balance constraint, node at battery, battery leakage modeled, battery efficiency modeled
        @constraint(m, [t=1:NumTime-1], InStorageBattery[t+1] == InStorageBattery[t] * δt * (1 - BatteryLoss) +  δt * (PV2B[t] * η - B2H[t])); # [kWh]

        # Battery discharging constraint, node at battery 
        @constraint(m, [t=1:NumTime], δt * B2H[t] <= InStorageBattery[t]); # [kWh]
        
        # Battery power inverter constraint, node at battery (inverter power constraint) 
        @constraint(m, [t=1:NumTime], B2H[t] + PV2B[t] <= InverterSize); # [kW]
        
        # Battery storage size constraint, node at battery
        @constraint(m, [t=1:NumTime], InStorageBattery[t] <= BatterySize); # [kWh]
        
        # Battery storage max discharge constraint, node at battery, always at least 20% full 
        @constraint(m, [t=1:NumTime], InStorageBattery[t] >= BatterySize * (1-MaxDischarge)); # [kWh]
       
        # Heating and Cooling band temperature controls
        begin 
            # Basically, all these constraints with all these binary variables are made to achieve one thing:
            #=
                The internal logic contains 3 data points “zone_temp_cooling_setpoint”, “zone_temp_heating_setpoint”, and “zone_temp_setpoint_delta”. 
                Based on the actual room temperature, these 3 data points determine whether the system should be in “heating”, “cooling” or “deadband” mode. 
                Let’s say these 3 data points are set to values of 273.15+22 (which is 22 °C), 273.15+20 (which is 20 °C), and 0.5, respectively. If the system 
                is originally in “heating” mode because the room temperature is 18 °C, then if the room temperature rises above the “zone_temp_heating_setpoint” 
                plus “zone_temp_setpoint_delta” (which is 20.5 °C), then the system will switch to “deadband” mode. When room temperature stays above 20.5 °C in “deadband” mode, 
                if the room temperature drops below “zone_temp_heating_setpoint” minus “zone_temp_setpoint_delta” (which is 19.5 °C), then the system will re-enter the “heating”. 

                This logic is the same for the “cooling” mode. If room temperature drops below 21.5 °C (“zone_temp_cooling_setpoint” minus “zone_temp_setpoint_delta”), 
                then the system will switch to “deadband” mode. If room temperature rises above 22.5 °C (“zone_temp_cooling_setpoint” plus “zone_temp_setpoint_delta”), 
                then the system will switch to “cooling” mode.  
            =#

            
            # Heating Logic
            begin
                #=
                # These are nonlinear (conditional) constraints that Gurobi cannot solve

                @constraint(m, [t=1:NumTime], heating_start[t] <= (TemperatureIndoor[t] <= zone_temp_heating_setpoint - zone_temp_setpoint_delta ? 1 : 0) * M); # [Bin]
                @constraint(m, [t=1:NumTime], heating_end[t] <= (TemperatureIndoor[t] >= zone_temp_heating_setpoint + zone_temp_setpoint_delta ? 1 : 0) * M); # [Bin]
                @constraint(m, [t=1:NumTime], heating_start[t] >= (TemperatureIndoor[t] <= zone_temp_heating_setpoint - zone_temp_setpoint_delta ? 1 : 0)); # [Bin]
                @constraint(m, [t=1:NumTime], heating_end[t] >= (TemperatureIndoor[t] >= zone_temp_heating_setpoint + zone_temp_setpoint_delta ? 1 : 0)); # [Bin]
                =#

                # These are the translated linear constraints
                @constraint(m, [t=1:NumTime], heating_start[t] <= heating_start_binary[t] * M)
                @constraint(m, [t=1:NumTime], heating_end[t] <= heating_end_binary[t] * M)
                @constraint(m, [t=1:NumTime], heating_start[t] >= heating_start_binary[t])
                @constraint(m, [t=1:NumTime], heating_end[t] >= heating_end_binary[t])

                # Relate binary variables to the temperature conditions using linear constraints
                @constraint(m, [t=1:NumTime], TemperatureIndoor[t] - (zone_temp_heating_setpoint - zone_temp_setpoint_delta) <= (1 - heating_start_binary[t]) * M)
                @constraint(m, [t=1:NumTime], (zone_temp_heating_setpoint - zone_temp_setpoint_delta) - TemperatureIndoor[t] <= heating_start_binary[t] * M)
                @constraint(m, [t=1:NumTime], TemperatureIndoor[t] - (zone_temp_heating_setpoint + zone_temp_setpoint_delta) >= (1 - heating_end_binary[t]) * -M)
                @constraint(m, [t=1:NumTime], (zone_temp_heating_setpoint + zone_temp_setpoint_delta) - TemperatureIndoor[t] >= heating_end_binary[t] * -M)

                # Ensure that heating[t] is 1 if heating_start[t] == 1 and heating_end[t] == 0
                @constraint(m, [t=1:NumTime], heating[t] >= heating_start[t] - heating_end[t]); # [Bin]
                
                # Additional constraints to ensure the logic is maintained over time
                @constraint(m, heating[1] <= heating_0 + heating_start[1]); # [Bin]
                @constraint(m, heating[1] <= 1 - heating_end[1]); # [Bin]
                @constraint(m, heating[1] >= heating_0 - heating_end[1]); # [Bin]

                @constraint(m, [t=2:NumTime], heating[t] <= heating[t-1] + heating_start[t]); # [Bin]
                @constraint(m, [t=2:NumTime], heating[t] <= 1 - heating_end[t]); # [Bin]
                @constraint(m, [t=2:NumTime], heating[t] >= heating[t-1] - heating_end[t]); # [Bin]
            end
            
            # Cooling Logic
            begin
                #=
                # These are nonlinear (conditional) constraints that Gurobi cannot solve

                @constraint(m, [t=1:NumTime], cooling_start[t] <= (TemperatureIndoor[t] >= zone_temp_cooling_setpoint + zone_temp_setpoint_delta ? 1 : 0) * M); # [Bin]
                @constraint(m, [t=1:NumTime], cooling_end[t] <= (TemperatureIndoor[t] <= zone_temp_cooling_setpoint - zone_temp_setpoint_delta ? 1 : 0) * M); # [Bin]
                @constraint(m, [t=1:NumTime], cooling_start[t] >= (TemperatureIndoor[t] >= zone_temp_cooling_setpoint + zone_temp_setpoint_delta ? 1 : 0)); # [Bin]
                @constraint(m, [t=1:NumTime], cooling_end[t] >= (TemperatureIndoor[t] <= zone_temp_cooling_setpoint - zone_temp_setpoint_delta ? 1 : 0)); # [Bin]
                =#

                # These are the translated linear constraints
                @constraint(m, [t=1:NumTime], cooling_start[t] <= cooling_start_binary[t] * M)
                @constraint(m, [t=1:NumTime], cooling_end[t] <= cooling_end_binary[t] * M)
                @constraint(m, [t=1:NumTime], cooling_start[t] >= cooling_start_binary[t])
                @constraint(m, [t=1:NumTime], cooling_end[t] >= cooling_end_binary[t])

                # Relate binary variables to the temperature conditions using linear constraints
                @constraint(m, [t=1:NumTime], TemperatureIndoor[t] - (zone_temp_cooling_setpoint + zone_temp_setpoint_delta) <= (1 - cooling_start_binary[t]) * M)
                @constraint(m, [t=1:NumTime], (zone_temp_cooling_setpoint + zone_temp_setpoint_delta) - TemperatureIndoor[t] <= cooling_start_binary[t] * M)
                @constraint(m, [t=1:NumTime], TemperatureIndoor[t] - (zone_temp_cooling_setpoint - zone_temp_setpoint_delta) >= (1 - cooling_end_binary[t]) * -M)
                @constraint(m, [t=1:NumTime], (zone_temp_cooling_setpoint - zone_temp_setpoint_delta) - TemperatureIndoor[t] >= cooling_end_binary[t] * -M)

                # Ensure that cooling[t] is 1 if cooling_start[t] == 1 and cooling_end[t] == 0
                @constraint(m, [t=1:NumTime], cooling[t] >= cooling_start[t] - cooling_end[t]); # [Bin]
                
                # Additional constraints to ensure the logic is maintained over time
                @constraint(m, cooling[1] <= cooling_0 + cooling_start[1]); # [Bin]
                @constraint(m, cooling[1] <= 1 - cooling_end[1]); # [Bin]
                @constraint(m, cooling[1] >= cooling_0 - cooling_end[1]); # [Bin]

                @constraint(m, [t=2:NumTime], cooling[t] <= cooling[t-1] + cooling_start[t]); # [Bin]
                @constraint(m, [t=2:NumTime], cooling[t] <= 1 - cooling_end[t]); # [Bin]
                @constraint(m, [t=2:NumTime], cooling[t] >= cooling[t-1] - cooling_end[t]); # [Bin]
            end
            
            @constraint(m, [t=1:NumTime], HP_state[2, t] + PCM_state[4, t] <= 2 * heating[t]); # [Bin]
            @constraint(m, [t=1:NumTime], HP_state[3, t] + PCM_state[5, t] <= 2 * cooling[t]); # [Bin]
            
            @constraint(m, [t=1:NumTime], HP_state[2, t] + PCM_state[4, t] >= heating[t]); # [Bin]
            @constraint(m, [t=1:NumTime], HP_state[3, t] + PCM_state[5, t] >= cooling[t]); # [Bin]

            # Heat Pump thermal energy balance constraints, node at Heat Pump
            @constraint(m, [t=1:NumTime], H2HP[t] * COP_H * 3412.14 - PCM_H_Charge[t] - HP2R[t] == 0); # [BTU/hr]
            @constraint(m, [t=1:NumTime], H2C[t] * COP_C * 3412.14 - PCM_C_Charge[t] - C2R[t] == 0); # [BTU/hr]
        end

        # PCM heating storage initialization constraint, node at PCM heating storage
        @constraint(m, PCM_SOC[1, 1] == InStoragePCM_H_1); # [1]

        # PCM cooling storage initialization constraint, node at PCM cooling storage
        @constraint(m, PCM_SOC[2, 1] == InStoragePCM_C_1); # [1]

        # PCM heating storage balance constraint, node at PCM heating storage
        @constraint(m, [t=1:NumTime-1], PCM_SOC[1, t+1] == PCM_SOC[1, t] + (δt * (PCM_H_Charge[t] - PCM_H_Discharge[t]))/PCM_H_Size); # [1]
        
        # PCM cooling storage balance constraint, node at PCM cooling storage
        @constraint(m, [t=1:NumTime-1], PCM_SOC[2, t+1] == PCM_SOC[2, t] + (δt * (PCM_C_Charge[t] - PCM_C_Discharge[t]))/PCM_C_Size); # [1]
        
        # PCM heating and cooling storage state of charge constraint, node at PCM storage
        @constraint(m, [ω = 1:2, t=2:NumTime], PCM_SOC[ω, t] <= 1.05); # [1]
    end
    ########### Solve  ##########
    optimize!(m);
    # Check the solver status
    status = JuMP.termination_status(m)
    println(status)
    println("Finished solving big model") 
    ########### Model Results  ##########
    begin
        # Collect J_States (the states that come from Julia MPC which will be used as input for the Julia MPC next iteration)
        J_States_new = zeros(3)
        J_States_new[1] = PVGeneration[1] * PVSize # [kW] PV power generation from time step 1
        J_States_new[2] = value.(E_total[1]) # [kW] Total Plug Load from time step 1
        J_States_new[3] = value.(InStorageBattery[1]) # [kWh] The battery SOC at the beginning of time step 1

        # Collect operational commands needed as input for FMU
        Modelica_Input = 0
        if value.(PCM_state[2, 1]) == 1
            Modelica_Input = 1
        elseif value.(PCM_state[3, 1]) == 1
            Modelica_Input = 2
        elseif value.(PCM_state[4, 1]) == 1 || value.(PCM_state[5, 1]) == 1
            Modelica_Input = 3
        else
            Modelica_Input = 0
        end
    end

    # Total penalty for loss of load
    TotalCost = objective_value(m) # [$]
    return TotalCost, J_States_new, [Modelica_Input, Modelica_Input]
end


#=
function Optimize1(iteration, Weather, Schedules, M_States, J_States)
    ########## Instructions  ##########
    # current optimization model output for modelica (0 to 3)
    # M_States should be a dataframe, which is the entire output of the FMU from previous iteration.

    ########## Data Preparations  ##########  
    begin
        # Set timesteps 
        NumTime = size(Weather)[1] # 168 hrs for now

        # Add forecast errors to the ambient temperature forecast
        Weather.Temperature = Add_Uncertainty_Temp(Weather.Temperature);

        Weather.CFM = zeros(NumTime)
        Weather.RadHeatGain = zeros(NumTime)
        Weather.RadCool = zeros(NumTime)
        Weather.PV = zeros(NumTime)

        for i = 1:NumTime
            Weather.CFM[i] = Al*sqrt(Cs*(T_indoor_constant + 273.15 - T_d) + Cw*wind_speed[i]^2) * 2.11888 # [ft^3/min] 
            Weather.RadHeatGain[i] = calculate_solarheatgain(Weather[i,:datetime], Weather[i,:DNI], Weather[i,:DHI], Weather[i,:GHI]) # [BTU/hr]
            day_status = IsDay(Weather.SolarTime[i], Weather.datetime[i]) 
            Weather.RadCool[i] = Q_radiativeCooling(T_indoor_constant, Weather.Temperature[i], Weather[i,:"Dew Point"], Weather.SolarTime[i], day_status)
            Weather.PV[i] = calculate_PV(Weather[i,:datetime], Weather[i,:DNI], Weather[i,:DHI], Weather[i,:GHI], 180-Berg_tilt, module_surface_tilt, Weather[i,:celltemp])
        end

        # Declare time series variables
        TemperatureAmbient = Weather.Temperature;
        PercentLighting = Schedules[:,1];
        PercentPlug = Schedules[:,2];
        PercentOccupied = Schedules[:,3];
        PVGeneration = Weather.PV;
        RadHeatGain = Weather.RadHeatGain;
        RadCooling = Weather.RadCool;
        CFM = Weather.CFM;

        # Length of the M_States dataframe
        msize = size(M_States)[1];

        InStoragePCM_C_1 = min(M_States[msize, :"PCM_Cold_SOC"], 1) # [1]
        InStoragePCM_H_1 = min(M_States[msize, :"PCM_Hot_SOC"], 1) # [1]
        TemperatureIndoor_1 = (M_States[msize, :"Indoor_Temp"] - 273.15) * (9/5) + 32 # [°F]
        InStorageBattery_1 = J_States[3] # [kWh]
        
        if iteration > 1
            println("M_States")
            print(M_States)
            println()
            time_intervals = zeros(msize-1); # [hr]

            for i = 1:msize-1
                time_intervals[i] = (M_States[i+1, :"Time"] - M_States[i, :"Time"])/3600 # [hr]
            end
            println("time intervals")
            print(time_intervals)
            println()
            println("E_Pump_1_0")
            E_Pump_1_0 = sum(((M_States[i+1, :"Pump 1 Electric Power"] + M_States[i, :"Pump 1 Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            print(E_Pump_1_0)
            println()
            E_Pump_2_0 = sum(((M_States[i+1, :"Pump 2 Electric Power"] + M_States[i, :"Pump 2 Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            E_HP_0 = sum(((M_States[i+1, :"HP Electric Power"] + M_States[i, :"HP Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            Thermal_Power_Delivered_0 = sum(((M_States[i+1, :"HP Useful Thermal Power"] + M_States[i, :"HP Useful Thermal Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            E_Fan_0 = sum(((M_States[i+1, :"Fan Coil Fan Electric Power"] + M_States[i, :"Fan Coil Fan Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            
            # Julia States
            PVGen_0 = J_States[1] # PV Generation at last timestep (Actual)
            E_other_0 = J_States[2] # [kW] Other Electrical Load (from schedule) at last timestep (Actual)
            B_SOC_0 = J_States[3] # Battery SOC at last timestep (Actual)
            
            E_Modelica_0 = E_Pump_1_0 + E_Pump_2_0 + E_HP_0 + E_Fan_0 # [kWh] Actual Electrical Load from HP system at last timestep
            Total_Load_0 = E_Modelica_0 + E_other_0 # [kWh] Actual total Electrical Load at last timestep
            
            # println("E_Modelica")
            # print(E_Modelica_0)
            # println()
            
            # Small Optimization Model to determine initial battery SOC
            begin
                println("Start solving small model")

                # Does the actual cost or loss of load come from here? Probably
                m_0 = Model(Gurobi.Optimizer)
                
                # define unit
                @variable(m_0, InStorageBattery_past[1:2] >= 0); # [kWh]
                @variable(m_0, Curtailment_past >= 0); # [kW]
                @variable(m_0, PV2B_past >= 0); # [kW]
                @variable(m_0, PV2H_past >= 0); # [kW]
                @variable(m_0, B2H_past >= 0); # [kW]
                @objective(m_0, Min, Curtailment_past); # [kW]

                @constraint(m_0, PVGen_0 == PV2B_past + PV2H_past + Curtailment_past) # [kW]
                @constraint(m_0, InStorageBattery_past[2] == InStorageBattery_past[1] * δt * (1 - BatteryLoss) +  δt * (PV2B_past * η - B2H_past)); # [kWh]
                @constraint(m_0, Total_Load_0 == PV2H_past * η_PVIV * δt + B2H_past * η * δt) # [kWh]

                # Battery discharging constraint, node at battery 
                @constraint(m_0, InStorageBattery_past[1] == B_SOC_0); # [kWh]

                # Battery discharging constraint, node at battery 
                @constraint(m_0, δt * B2H_past <= InStorageBattery_past[1]); # [kWh]
                
                # Battery power inverter constraint, node at battery (inverter power constraint) 
                @constraint(m_0, B2H_past + PV2B_past <= InverterSize); # [kW]
                
                # Battery storage size constraint, node at battery
                @constraint(m_0, [t=1:2], InStorageBattery_past[t] <= BatterySize); # [kWh]
                
                # Battery storage max discharge constraint, node at battery, always at least 20% full 
                @constraint(m_0, [t=1:2], InStorageBattery_past[t] >= BatterySize * (1-MaxDischarge)); # [kWh]

                optimize!(m_0);
                
                # Actual Battery SOC to start with
                InStorageBattery_1 = value.(InStorageBattery_past[2]) # [kWh]
                println("Finished solving small model")
                #=
                println()
                println("Battery SOC is:")
                print(InStorageBattery_1)
                println()
                println("Hot PCM SOC is:")
                print(InStoragePCM_H_1)
                println()
                println("Cold PCM SOC is:")
                print(InStoragePCM_C_1)
                println()
                println("Indoor Temperature is:")
                print(TemperatureIndoor_1)
                println()
                =#
            end
        end    
    end
    println("Start solving big model")  
    ########## Declare model  ##########
    begin
        # Define the model name and solver. In this case, model name is "m"
        # m = Model(Clp.Optimizer)
        # m = Model(Ipopt.Optimizer)
        begin
            # Set path to license (for those using Gurobi)
            # ENV["GRB_LICENSE_FILE"] = "/Library/gurobi1100/gurobi.lic"
            m = Model(Gurobi.Optimizer)
        end
    end
    ######## Decision variables ########
    begin
        @variable(m, PV2H[1:NumTime] >= 0); # [kW] electrical power transfer from PV to home (Berg)

        @variable(m, PV2G[1:NumTime] >= 0); # [kW] electrical power transfer from PV to ground (curtailment)

        @variable(m, PV2B[1:NumTime] >= 0); # [kW] electrical power transfer from PV to battery

        @variable(m, B2H[1:NumTime] >= 0); # [kW] electrical power transfer from battery to home (Berg)

        @variable(m, H2HP[1:NumTime] >= 0); # [kW] electrical power transfer from home (Berg) to heat pump heating mode

        @variable(m, HP2H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from heat pump heating mode to home (Berg)

        @variable(m, H2C[1:NumTime] >= 0); # [kW] electrical power transfer from home (Berg) to heat pump cooling unit

        @variable(m, C2H[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from heat pump cooling mode to home (Berg)
        
        @variable(m, HP_OPTION[1:NumTime], Bin) # Heating Mode = 1, Cooling Mode = 0

        @variable(m, TES_Charging_OPTION[1:NumTime], Bin) # Charge hot TES = 1, Charge cold TES = 0

        @variable(m, TES_Discharging_OPTION[1:NumTime], Bin) # Discharge hot TES = 1, Discharge cold TES = 0

        @variable(m, TES_CD_OPTION[1:NumTime], Bin) # Charge TES = 1, Discharge TES = 0

        @variable(m, B_OPTION[1:NumTime], Bin) # Charge battery = 1, Discharge battery = 0

        @variable(m, HP2PCM_H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from heat pump heating mode to PCM heating storage

        @variable(m, C2PCM_C[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from heat pump cooling mode to PCM cooling storage

        @variable(m, PCM_H2H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from PCM heating storage to home (Berg)

        @variable(m, PCM_C2H[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from PCM cooling storage to home (Berg)

        @variable(m, InStorageBattery[1:NumTime] >= 0); # [kWh] Battery Remaining Charge

        @variable(m, InStoragePCM_H[1:NumTime] >= 0); # [BTU] PCM Heating Remaining Charge

        @variable(m, InStoragePCM_C[1:NumTime] >= 0); # [BTU] PCM Cooling Remaining Charge

        @variable(m, TemperatureIndoor[1:NumTime] >= 0); # [°F] Indoor Air Temperature

        @variable(m, G2H[1:NumTime] >= 0); # [kW] Loss of Load
    end
    ############ Objective Functions #############
    begin
        # Set single objective for minimizing annual total cost
        
        # Calculate Yearly Short Run Marginal Cost [$]
        # @expression(m, short_run_marginal_cost, δt * sum(B2H[t] * C_B_OPV for t=1:NumTime))
        
        # Penalty for curtailment [$]
        # @expression(m, penalty_c, δt * sum(PV2G[t] for t = 1:NumTime))

        # Penalty for loss of load [$]
        @expression(m, penalty_l, M * δt * sum(G2H[t] for t = 1:NumTime))

        # Total Cost over Optimization Horizon [$]
        # @objective(m, Min, short_run_marginal_cost + penalty_c + penalty_l);
        @objective(m, Min, δt * sum(G2H[t] for t = 1:NumTime));
    end
    ############# Expressions ############
    begin
        # DeltaTemp 
        @expression(m, TempDelta[t=1:NumTime], TemperatureAmbient[t] - TemperatureIndoor[t]); # [°F]

        # DeltaTemp 
        @expression(m, TempDelta_C[t=1:NumTime], (TemperatureAmbient[t] * 1.8 + 32) - (TemperatureIndoor[t] * 1.8 + 32)); # [°C]

        # Electricity usage from lighting 
        @expression(m, E_Lighting[t=1:NumTime], PeakLighting * PercentLighting[t]); # [kW]

        # Electricity usage from plugs
        @expression(m, E_Plugs[t=1:NumTime], PeakPlugLoad * PercentPlug[t]); # [kW]

        # Total electricity usage
        @expression(m, E_total[t=1:NumTime], E_Lighting[t] + E_Plugs[t]); # [kW]

        # Calculate Ventilation
        @expression(m, CFMVen[t=1:NumTime], min(Rp * PercentOccupied[t] * MaxOccupancy + Ra * Area, PercentOccupied[t] * MaxOccupancy * Ventilation)); # [ft^3/min]

        # Heat gain from lighting
        @expression(m, Q_Lighting[t=1:NumTime], PeakLighting * PercentLighting[t] * 3412.14); # [BTU/hr]

        # Heat gain from plugs
        @expression(m, Q_Plugs[t=1:NumTime], PeakPlugLoad * PercentPlug[t] * 3412.14); # [BTU/hr]

        # Heat gain from occupancy
        @expression(m, Q_Occupancy[t=1:NumTime], MaxOccupancy * PercentOccupied[t] * TotalPersonHeat); # [BTU/hr]

        # Heat gain from infiltration
        # 1.08 = specific heat capacity of air at STP:0.24 [BTU/(lb*°F)] * air density at STP:0.075 [lb/ft^3] * 60 [min/hr] 
        @expression(m, Q_Infiltration[t=1:NumTime], 1.08 * TempDelta[t] * CFM[t]); # [BTU/hr]

        # Heat gain from ventilation
        @expression(m, Q_Ventilation[t=1:NumTime], 1.08 * TempDelta[t] * CFMVen[t]); # [BTU/hr]

        # Heat gain from lighting, plugs, occupancy, ventilation, infiltration
        @expression(m, Q_Others[t=1:NumTime], Q_Infiltration[t] + Q_Ventilation[t] + Q_Occupancy[t] + Q_Lighting[t] + Q_Plugs[t]); # [BTU/hr]

        # Heat gain through structural evenlope
        @expression(m, Q_Envelope[t=1:NumTime], UA * TempDelta[t]); # [BTU/hr]

        # Heat gain from solar radiation
        @expression(m, Q_Rad[t=1:NumTime], SHGC * RadHeatGain[t]); # [BTU/hr]

        # Radiative cooling
        @expression(m, Q_RadCool[t=1:NumTime], RCC * RadCooling[t]); # [BTU/hr]

        # Detailed COP of heating and cooling (requires nonlinear optimization solver, although we can hold indoor temperature as a constant to keep the model linear)
        # COP (Heating Home)
        # @expression(m, COP_Heating_H[t=1:NumTime], HP_a + HP_b * (-TempDelta_C[t]) + HP_c * TempDelta[t]^2);
        
        # COP (Heating PCM H)
        # @expression(m, COP_Heating_PCM[t=1:NumTime], HP_a + HP_b * (48 - Weather[t, 9]) + HP_c * (48 - Weather[t, 9])^2); 

        # COP (Cooling Home)
        # @expression(m, COP_Cooling_H[t=1:NumTime], HP_a + HP_b * (TempDelta_C[t]) + HP_c * TempDelta[t]^2);
        
        # COP (Cooling PCM C)
        # @expression(m, COP_Cooling_PCM[t=1:NumTime], HP_a + HP_b * (Weather[t, 9] - 11) + HP_c * (Weather[t, 9] - 11)^2); 
    end
    ############# Constraints ############
    begin
        # Initialization Constraints with States from last timestep
        begin
            # Set point temperature range constraints
            @constraint(m, TemperatureIndoor[1] == TemperatureIndoor_1); # [°F]

            # PCM heating storage initialization constraint, node at PCM heating storage
            @constraint(m, InStoragePCM_H[1] == InStoragePCM_H_1 * PCM_H_Size); # [BTU]

            # PCM cooling storage initialization constraint, node at PCM cooling storage
            @constraint(m, InStoragePCM_C[1] == InStoragePCM_C_1 * PCM_C_Size); # [BTU]

            # Battery storage initialization constraint
            @constraint(m, InStorageBattery[1] == InStorageBattery_1); # [kWh] 
        end

        # Set point temperature range constraints 
        @constraint(m, [t=1:NumTime], TemperatureIndoor[t] <= SetPointT_High); # [°F]

        @constraint(m, [t=1:NumTime], TemperatureIndoor[t] >= SetPointT_Low); # [°F]

        # Internal temperature balance evolution constraint
        @constraint(m, [t=1:NumTime-1], TemperatureIndoor[t+1] == TemperatureIndoor[t] + 
        (δt/TC)*(Q_Others[t] + Q_Envelope[t] + Q_Rad[t] - Q_RadCool[t] + HP2H[t] - C2H[t] + PCM_H2H[t] - PCM_C2H[t])); # [°F]
        
        # PV energy balance constraint, node at PV
        @constraint(m, [t=1:NumTime], PVGeneration[t] * PVSize ==  PV2B[t] + PV2H[t] + PV2G[t]); # [kW]
    
        # House electricity load constraint, node at house, battery efficiency modeled
        @constraint(m, [t=1:NumTime], E_total[t] + H2HP[t] + H2C[t] == PV2H[t] * η_PVIV + B2H[t] * η + G2H[t]); # [kW]

        # Battery storage balance constraint, node at battery, battery leakage modeled, battery efficiency modeled
        @constraint(m, [t=1:NumTime-1], InStorageBattery[t+1] == InStorageBattery[t] * δt * (1 - BatteryLoss) +  δt * (PV2B[t] * η - B2H[t])); # [kWh]

        # Battery discharging constraint, node at battery 
        @constraint(m, [t=1:NumTime], δt * B2H[t] <= InStorageBattery[t]); # [kWh]
        
        # Battery power inverter constraint, node at battery (inverter power constraint) 
        @constraint(m, [t=1:NumTime], B2H[t] + PV2B[t] <= InverterSize); # [kW]
        
        # Battery storage size constraint, node at battery
        @constraint(m, [t=1:NumTime], InStorageBattery[t] <= BatterySize); # [kWh]
        
        # Battery storage max discharge constraint, node at battery, always at least 20% full 
        @constraint(m, [t=1:NumTime], InStorageBattery[t] >= BatterySize * (1-MaxDischarge)); # [kWh]
        
        # Heating and Cooling Constraints (Needs Remodeling)
        begin
            #=
            # Detailed COP of heating and cooling (requires nonlinear optimization solver)
            # Heating (from kWh to BTU), node at heat pump heating option
            @constraint(m, [t=1:NumTime], H2HP[t] * 3412.14 == HP2PCM_H[t]/COP_Heating_PCM[t] + HP2H[t]/COP_Heating_H[t]);
            
            # Cooling (from kWh to BTU), node at heat pump cooling option
            @constraint(m, [t=1:NumTime], H2C[t] * 3412.14 == C2PCM_C[t]/COP_Cooling_PCM[t] + C2H[t]/COP_Cooling_H[t]);
            =#

            # Heating (from kW to BTU/hr), node at heat pump heating option
            @constraint(m, [t=1:NumTime], H2HP[t] * 3412.14 == (HP2PCM_H[t] + HP2H[t])/COP_H); # [BTU/hr]

            # Cooling (from kW to BTU/hr), node at heat pump cooling option
            @constraint(m, [t=1:NumTime], H2C[t] * 3412.14 == (C2PCM_C[t] + C2H[t])/COP_C); # [BTU/hr]

            # Heat pump power capacity constraint, node at heat pump heating mode
            @constraint(m, [t=1:NumTime], H2HP[t] <= Cap_HP_H); # [kW]

            # Heat pump power capacity constraint, node at heat pump cooling mode
            @constraint(m, [t=1:NumTime], H2C[t] <= Cap_HP_C); # [kW]

            # PCM heating storage balance constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime-1], InStoragePCM_H[t+1] == InStoragePCM_H[t] + δt * (HP2PCM_H[t] - PCM_H2H[t])); # [BTU]
            
            # PCM cooling storage balance constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime-1], InStoragePCM_C[t+1] == InStoragePCM_C[t] + δt * (C2PCM_C[t] - PCM_C2H[t])); # [BTU]
            
            # PCM heating storage discharging constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime], δt * PCM_H2H[t] <= InStoragePCM_H[t]); # [BTU]
            
            # PCM cooling storage discharging constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime], δt * PCM_C2H[t] <= InStoragePCM_C[t]); # [BTU]
            
            # PCM heating storage size constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime], InStoragePCM_H[t] <= PCM_H_Size); # [BTU]

            # PCM cooling storage size constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime], InStoragePCM_C[t] <= PCM_C_Size); # [BTU]
        end
        
        # Single Operation Choice Constraints
        begin
            # If HP_OPTION = 1, heating mode on, heating only
            @constraint(m, [t=1:NumTime], H2HP[t] <= M * HP_OPTION[t]) # [kW]          

            # If HP_OPTION = 0, cooling mode on, cooling only
            @constraint(m, [t=1:NumTime], H2C[t] <= M * (1 - HP_OPTION[t])) # [kW]    
    
            # If TES_Charging_OPTION = 1, charge hot TES only
            @constraint(m, [t=1:NumTime], HP2PCM_H[t] <= M * TES_Charging_OPTION[t]) # [BTU/hr]     
    
            # If TES_Charging_OPTION = 0, charge cold TES only
            @constraint(m, [t=1:NumTime], C2PCM_C[t] <= M * (1 - TES_Charging_OPTION[t])) # [BTU/hr]

            # If TES_Discharging_OPTION = 1, discharge hot TES only
            @constraint(m, [t=1:NumTime], PCM_H2H[t] <= M * TES_Discharging_OPTION[t]) # [BTU/hr]          
    
            # If TES_Discharging_OPTION = 0, discharge cold TES only
            @constraint(m, [t=1:NumTime], PCM_C2H[t] <= M * (1 - TES_Discharging_OPTION[t])) # [BTU/hr]

            # Total TES Charging 
            @expression(m, TES_Charge[t=1:NumTime], HP2PCM_H[t] + C2PCM_C[t]); # [BTU/hr]

            # Total TES Discharging 
            @expression(m, TES_Discharge[t=1:NumTime], PCM_H2H[t] + PCM_C2H[t]); # [BTU/hr]

            # If TES_CD_OPTION = 1, charge only
            @constraint(m, [t=1:NumTime], TES_Charge[t] <= M * TES_CD_OPTION[t]) # [BTU/hr]          

            # If TES_CD_OPTION = 0, discharge only
            @constraint(m, [t=1:NumTime], TES_Discharge[t] <= M * (1 - TES_CD_OPTION[t])) # [BTU/hr] 

            # If B_OPTION = 1, charge only
            @constraint(m, [t=1:NumTime], PV2B[t] <= M * B_OPTION[t]) # [kW]          

            # If B_OPTION = 0, discharge only
            @constraint(m, [t=1:NumTime], B2H[t] <= M * (1 - B_OPTION[t])) # [kW] 
        end

    end
    ########### Solve  ##########
    optimize!(m);
    println("Finished solving big model") 
    ########### Model Results  ##########
    begin
        # Collect J_States (the states that come from Julia MPC which will be used as input for the Julia MPC next iteration)
        J_States_new = zeros(3)
        J_States_new[1] = PVGeneration[1] * PVSize # [kW] PV power generation from time step 1
        J_States_new[2] = value.(E_total[1]) # [kW] Total Plug Load from time step 1
        J_States_new[3] = value.(InStorageBattery[1]) # [kWh] The battery SOC at the beginning of time step 1

        # Collect operational commands needed as input for FMU
        Modelica_Input = 0
        if value.(TES_Discharge[1]) > 0
            Modelica_Input = 3
        elseif value.(HP2PCM_H[1]) > 0
            Modelica_Input = 1
        elseif value.(C2PCM_C[1]) > 0
            Modelica_Input = 2
        end
    end
    
    # Total penalty for loss of load
    TotalCost = objective_value(m) # [$]
    return TotalCost, J_States_new, [Modelica_Input, Modelica_Input]
end


function Optimize3(iteration, Weather, Schedules, M_States, J_States)
    ########## Instructions  ##########
    # current optimization model output for modelica (0 to 3)
    # M_States should be a dataframe, which is the entire output of the FMU from previous iteration.

    ########## Data Preparations  ##########  
    begin
        # Set timesteps 
        NumTime = size(Weather)[1] # 168 hrs for now

        # Add forecast errors to the ambient temperature forecast
        Weather.Temperature = Add_Uncertainty_Temp(Weather.Temperature);

        Weather.CFM = zeros(NumTime)
        Weather.RadHeatGain = zeros(NumTime)
        Weather.RadCool = zeros(NumTime)
        Weather.PV = zeros(NumTime)

        for i = 1:NumTime
            Weather.CFM[i] = Al*sqrt(Cs*(T_indoor_constant + 273.15 - T_d) + Cw*wind_speed[i]^2) * 2.11888 # [ft^3/min] 
            Weather.RadHeatGain[i] = calculate_solarheatgain(Weather[i,:datetime], Weather[i,:DNI], Weather[i,:DHI], Weather[i,:GHI]) # [BTU/hr]
            day_status = IsDay(Weather.SolarTime[i], Weather.datetime[i]) 
            Weather.RadCool[i] = Q_radiativeCooling(T_indoor_constant, Weather.Temperature[i], Weather[i,:"Dew Point"], Weather.SolarTime[i], day_status)
            Weather.PV[i] = calculate_PV(Weather[i,:datetime], Weather[i,:DNI], Weather[i,:DHI], Weather[i,:GHI], 180-Berg_tilt, module_surface_tilt, Weather[i,:celltemp])
        end

        # Declare time series variables
        TemperatureAmbient = Weather.Temperature;
        PercentLighting = Schedules[:,1];
        PercentPlug = Schedules[:,2];
        PercentOccupied = Schedules[:,3];
        PVGeneration = Weather.PV;
        RadHeatGain = Weather.RadHeatGain;
        RadCooling = Weather.RadCool;
        CFM = Weather.CFM;

        # Length of the M_States dataframe
        msize = size(M_States)[1];

        InStoragePCM_C_1 = min(M_States[msize, :"PCM_Cold_SOC"], 1) # [1]
        InStoragePCM_H_1 = min(M_States[msize, :"PCM_Hot_SOC"], 1) # [1]
        InStorageBattery_1 = J_States[3] # [kWh]
        TemperatureIndoor_1 = (M_States[msize, :"Indoor_Temp"] - 273.15) * (9/5) + 32 # [°F]

        # Variation setpoint for temperature
        begin
            # The time allowed to have variable set point temperatures before returning back to the original set point temperature bounds
            Temp_Flexible_Time = 2 # [hr]

            SetPointT_High_Variable = SetPointT_High
            SetPointT_Low_Variable = SetPointT_Low
            if TemperatureIndoor_1 >= SetPointT_High
                SetPointT_High_Variable = TemperatureIndoor_1
            end
    
            if TemperatureIndoor_1 <= SetPointT_Low
                SetPointT_Low_Variable = TemperatureIndoor_1
            end
        end

        
        if iteration > 1
            println("M_States")
            print(M_States)
            println()
            time_intervals = zeros(msize-1); # [hr]

            for i = 1:msize-1
                time_intervals[i] = (M_States[i+1, :"Time"] - M_States[i, :"Time"])/3600 # [hr]
            end
            println("time intervals")
            print(time_intervals)
            println()
            println("E_Pump_1_0")
            E_Pump_1_0 = sum(((M_States[i+1, :"Pump 1 Electric Power"] + M_States[i, :"Pump 1 Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            print(E_Pump_1_0)
            println()
            E_Pump_2_0 = sum(((M_States[i+1, :"Pump 2 Electric Power"] + M_States[i, :"Pump 2 Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            E_HP_0 = sum(((M_States[i+1, :"HP Electric Power"] + M_States[i, :"HP Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            Thermal_Power_Delivered_0 = sum(((M_States[i+1, :"HP Useful Thermal Power"] + M_States[i, :"HP Useful Thermal Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            E_Fan_0 = sum(((M_States[i+1, :"Fan Coil Fan Electric Power"] + M_States[i, :"Fan Coil Fan Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
            
            # Julia States
            PVGen_0 = J_States[1] # PV Generation at last timestep (Actual)
            E_other_0 = J_States[2] # [kW] Other Electrical Load (from schedule) at last timestep (Actual)
            B_SOC_0 = J_States[3] # Battery SOC at last timestep (Actual)
            
            E_Modelica_0 = E_Pump_1_0 + E_Pump_2_0 + E_HP_0 + E_Fan_0 # [kWh] Actual Electrical Load from HP system at last timestep
            Total_Load_0 = E_Modelica_0 + E_other_0 # [kWh] Actual total Electrical Load at last timestep
            
            # println("E_Modelica")
            # print(E_Modelica_0)
            # println()
            
            # Small Optimization Model to determine initial battery SOC
            begin
                println("Start solving small model")

                # Does the actual cost or loss of load come from here? Probably
                m_0 = Model(Gurobi.Optimizer)
                
                # define unit
                @variable(m_0, InStorageBattery_past[1:2] >= 0); # [kWh]
                @variable(m_0, Curtailment_past >= 0); # [kW]
                @variable(m_0, PV2B_past >= 0); # [kW]
                @variable(m_0, PV2H_past >= 0); # [kW]
                @variable(m_0, B2H_past >= 0); # [kW]
                @objective(m_0, Min, Curtailment_past); # [kW]

                @constraint(m_0, PVGen_0 == PV2B_past + PV2H_past + Curtailment_past) # [kW]
                @constraint(m_0, InStorageBattery_past[2] == InStorageBattery_past[1] * δt * (1 - BatteryLoss) +  δt * (PV2B_past * η - B2H_past)); # [kWh]
                @constraint(m_0, Total_Load_0 == PV2H_past * η_PVIV * δt + B2H_past * η * δt) # [kWh]

                # Battery discharging constraint, node at battery 
                @constraint(m_0, InStorageBattery_past[1] == B_SOC_0); # [kWh]

                # Battery discharging constraint, node at battery 
                @constraint(m_0, δt * B2H_past <= InStorageBattery_past[1]); # [kWh]
                
                # Battery power inverter constraint, node at battery (inverter power constraint) 
                @constraint(m_0, B2H_past + PV2B_past <= InverterSize); # [kW]
                
                # Battery storage size constraint, node at battery
                @constraint(m_0, [t=1:2], InStorageBattery_past[t] <= BatterySize); # [kWh]
                
                # Battery storage max discharge constraint, node at battery, always at least 20% full 
                @constraint(m_0, [t=1:2], InStorageBattery_past[t] >= BatterySize * (1-MaxDischarge)); # [kWh]

                optimize!(m_0);
                
                # Actual Battery SOC to start with
                InStorageBattery_1 = value.(InStorageBattery_past[2]) # [kWh]
                println("Finished solving small model")
                #=
                println()
                println("Battery SOC is:")
                print(InStorageBattery_1)
                println()
                println("Hot PCM SOC is:")
                print(InStoragePCM_H_1)
                println()
                println("Cold PCM SOC is:")
                print(InStoragePCM_C_1)
                println()
                println("Indoor Temperature is:")
                print(TemperatureIndoor_1)
                println()
                =#
            end
        end    
    end
    println("Start solving big model")  
    ########## Declare model  ##########
    begin
        # Define the model name and solver. In this case, model name is "m"
        # m = Model(Clp.Optimizer)
        # m = Model(Ipopt.Optimizer)
        begin
            # Set path to license (for those using Gurobi)
            # ENV["GRB_LICENSE_FILE"] = "/Library/gurobi1100/gurobi.lic"
            m = Model(Gurobi.Optimizer)
        end
    end
    ######## Decision variables ########
    begin
        @variable(m, PV2H[1:NumTime] >= 0); # [kW] electrical power transfer from PV to home (Berg)

        @variable(m, PV2G[1:NumTime] >= 0); # [kW] electrical power transfer from PV to ground (curtailment)

        @variable(m, PV2B[1:NumTime] >= 0); # [kW] electrical power transfer from PV to battery

        @variable(m, B2H[1:NumTime] >= 0); # [kW] electrical power transfer from battery to home (Berg)

        @variable(m, H2HP[1:NumTime] >= 0); # [kW] electrical power transfer from home (Berg) to heat pump heating mode

        @variable(m, HP2H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from heat pump heating mode to home (Berg)

        @variable(m, H2C[1:NumTime] >= 0); # [kW] electrical power transfer from home (Berg) to heat pump cooling unit

        @variable(m, C2H[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from heat pump cooling mode to home (Berg)
        
        @variable(m, HP_OPTION[1:NumTime], Bin) # Heating Mode = 1, Cooling Mode = 0

        @variable(m, TES_Charging_OPTION[1:NumTime], Bin) # Charge hot TES = 1, Charge cold TES = 0

        @variable(m, TES_Discharging_OPTION[1:NumTime], Bin) # Discharge hot TES = 1, Discharge cold TES = 0

        @variable(m, TES_CD_OPTION[1:NumTime], Bin) # Charge TES = 1, Discharge TES = 0

        @variable(m, B_OPTION[1:NumTime], Bin) # Charge battery = 1, Discharge battery = 0

        @variable(m, HP2PCM_H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from heat pump heating mode to PCM heating storage

        @variable(m, C2PCM_C[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from heat pump cooling mode to PCM cooling storage

        @variable(m, PCM_H2H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from PCM heating storage to home (Berg)

        @variable(m, PCM_C2H[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from PCM cooling storage to home (Berg)

        @variable(m, InStorageBattery[1:NumTime] >= 0); # [kWh] Battery Remaining Charge

        @variable(m, InStoragePCM_H[1:NumTime] >= 0); # [BTU] PCM Heating Remaining Charge

        @variable(m, InStoragePCM_C[1:NumTime] >= 0); # [BTU] PCM Cooling Remaining Charge

        @variable(m, TemperatureIndoor[1:NumTime] >= 0); # [°F] Indoor Air Temperature

        @variable(m, G2H[1:NumTime] >= 0); # [kW] Loss of Load
    end
    ############ Objective Functions #############
    begin
        # Set single objective for minimizing annual total cost
        
        # Calculate Yearly Short Run Marginal Cost [$]
        # @expression(m, short_run_marginal_cost, δt * sum(B2H[t] * C_B_OPV for t=1:NumTime))
        
        # Penalty for curtailment [$]
        # @expression(m, penalty_c, δt * sum(PV2G[t] for t = 1:NumTime))

        # Penalty for loss of load [$]
        @expression(m, penalty_l, M * δt * sum(G2H[t] for t = 1:NumTime))

        # Total Cost over Optimization Horizon [$]
        # @objective(m, Min, short_run_marginal_cost + penalty_c + penalty_l);
        @objective(m, Min, δt * sum(G2H[t] for t = 1:NumTime));
    end
    ############# Expressions ############
    begin
        # DeltaTemp 
        @expression(m, TempDelta[t=1:NumTime], TemperatureAmbient[t] - TemperatureIndoor[t]); # [°F]

        # DeltaTemp 
        @expression(m, TempDelta_C[t=1:NumTime], (TemperatureAmbient[t] * 1.8 + 32) - (TemperatureIndoor[t] * 1.8 + 32)); # [°C]

        # Electricity usage from lighting 
        @expression(m, E_Lighting[t=1:NumTime], PeakLighting * PercentLighting[t]); # [kW]

        # Electricity usage from plugs
        @expression(m, E_Plugs[t=1:NumTime], PeakPlugLoad * PercentPlug[t]); # [kW]

        # Total electricity usage
        @expression(m, E_total[t=1:NumTime], E_Lighting[t] + E_Plugs[t]); # [kW]

        # Calculate Ventilation
        @expression(m, CFMVen[t=1:NumTime], min(Rp * PercentOccupied[t] * MaxOccupancy + Ra * Area, PercentOccupied[t] * MaxOccupancy * Ventilation)); # [ft^3/min]

        # Heat gain from lighting
        @expression(m, Q_Lighting[t=1:NumTime], PeakLighting * PercentLighting[t] * 3412.14); # [BTU/hr]

        # Heat gain from plugs
        @expression(m, Q_Plugs[t=1:NumTime], PeakPlugLoad * PercentPlug[t] * 3412.14); # [BTU/hr]

        # Heat gain from occupancy
        @expression(m, Q_Occupancy[t=1:NumTime], MaxOccupancy * PercentOccupied[t] * TotalPersonHeat); # [BTU/hr]

        # Heat gain from infiltration
        # 1.08 = specific heat capacity of air at STP:0.24 [BTU/(lb*°F)] * air density at STP:0.075 [lb/ft^3] * 60 [min/hr] 
        @expression(m, Q_Infiltration[t=1:NumTime], 1.08 * TempDelta[t] * CFM[t]); # [BTU/hr]

        # Heat gain from ventilation
        @expression(m, Q_Ventilation[t=1:NumTime], 1.08 * TempDelta[t] * CFMVen[t]); # [BTU/hr]

        # Heat gain from lighting, plugs, occupancy, ventilation, infiltration
        @expression(m, Q_Others[t=1:NumTime], Q_Infiltration[t] + Q_Ventilation[t] + Q_Occupancy[t] + Q_Lighting[t] + Q_Plugs[t]); # [BTU/hr]

        # Heat gain through structural evenlope
        @expression(m, Q_Envelope[t=1:NumTime], UA * TempDelta[t]); # [BTU/hr]

        # Heat gain from solar radiation
        @expression(m, Q_Rad[t=1:NumTime], SHGC * RadHeatGain[t]); # [BTU/hr]

        # Radiative cooling
        @expression(m, Q_RadCool[t=1:NumTime], RCC * RadCooling[t]); # [BTU/hr]

        # Detailed COP of heating and cooling (requires nonlinear optimization solver, although we can hold indoor temperature as a constant to keep the model linear)
        # COP (Heating Home)
        # @expression(m, COP_Heating_H[t=1:NumTime], HP_a + HP_b * (-TempDelta_C[t]) + HP_c * TempDelta[t]^2);
        
        # COP (Heating PCM H)
        # @expression(m, COP_Heating_PCM[t=1:NumTime], HP_a + HP_b * (48 - Weather[t, 9]) + HP_c * (48 - Weather[t, 9])^2); 

        # COP (Cooling Home)
        # @expression(m, COP_Cooling_H[t=1:NumTime], HP_a + HP_b * (TempDelta_C[t]) + HP_c * TempDelta[t]^2);
        
        # COP (Cooling PCM C)
        # @expression(m, COP_Cooling_PCM[t=1:NumTime], HP_a + HP_b * (Weather[t, 9] - 11) + HP_c * (Weather[t, 9] - 11)^2); 
    end
    ############# Constraints ############
    begin
        # Initialization Constraints with States from last timestep
        begin
            # Set point temperature range constraints
            @constraint(m, TemperatureIndoor[1] == TemperatureIndoor_1); # [°F]

            # PCM heating storage initialization constraint, node at PCM heating storage
            @constraint(m, InStoragePCM_H[1] == InStoragePCM_H_1 * PCM_H_Size); # [BTU]

            # PCM cooling storage initialization constraint, node at PCM cooling storage
            @constraint(m, InStoragePCM_C[1] == InStoragePCM_C_1 * PCM_C_Size); # [BTU]

            # Battery storage initialization constraint
            @constraint(m, InStorageBattery[1] == InStorageBattery_1); # [kWh] 
        end

        # Set point temperature range constraints 
        @constraint(m, [t=1:NumTime], TemperatureIndoor[t] <= SetPointT_High); # [°F]

        @constraint(m, [t=1:NumTime], TemperatureIndoor[t] >= SetPointT_Low); # [°F]

        # Internal temperature balance evolution constraint
        @constraint(m, [t=1:NumTime-1], TemperatureIndoor[t+1] == TemperatureIndoor[t] + 
        (δt/TC)*(Q_Others[t] + Q_Envelope[t] + Q_Rad[t] - Q_RadCool[t] + HP2H[t] - C2H[t] + PCM_H2H[t] - PCM_C2H[t])); # [°F]
        
        # PV energy balance constraint, node at PV
        @constraint(m, [t=1:NumTime], PVGeneration[t] * PVSize ==  PV2B[t] + PV2H[t] + PV2G[t]); # [kW]
    
        # House electricity load constraint, node at house, battery efficiency modeled
        @constraint(m, [t=1:NumTime], E_total[t] + H2HP[t] + H2C[t] == PV2H[t] * η_PVIV + B2H[t] * η + G2H[t]); # [kW]

        # Battery storage balance constraint, node at battery, battery leakage modeled, battery efficiency modeled
        @constraint(m, [t=1:NumTime-1], InStorageBattery[t+1] == InStorageBattery[t] * δt * (1 - BatteryLoss) +  δt * (PV2B[t] * η - B2H[t])); # [kWh]

        # Battery discharging constraint, node at battery 
        @constraint(m, [t=1:NumTime], δt * B2H[t] <= InStorageBattery[t]); # [kWh]
        
        # Battery power inverter constraint, node at battery (inverter power constraint) 
        @constraint(m, [t=1:NumTime], B2H[t] + PV2B[t] <= InverterSize); # [kW]
        
        # Battery storage size constraint, node at battery
        @constraint(m, [t=1:NumTime], InStorageBattery[t] <= BatterySize); # [kWh]
        
        # Battery storage max discharge constraint, node at battery, always at least 20% full 
        @constraint(m, [t=1:NumTime], InStorageBattery[t] >= BatterySize * (1-MaxDischarge)); # [kWh]
        
        # Heating and Cooling Constraints (Needs Remodeling)
        begin
            #=
            # Detailed COP of heating and cooling (requires nonlinear optimization solver)
            # Heating (from kWh to BTU), node at heat pump heating option
            @constraint(m, [t=1:NumTime], H2HP[t] * 3412.14 == HP2PCM_H[t]/COP_Heating_PCM[t] + HP2H[t]/COP_Heating_H[t]);
            
            # Cooling (from kWh to BTU), node at heat pump cooling option
            @constraint(m, [t=1:NumTime], H2C[t] * 3412.14 == C2PCM_C[t]/COP_Cooling_PCM[t] + C2H[t]/COP_Cooling_H[t]);
            =#

            # Heating (from kW to BTU/hr), node at heat pump heating option
            @constraint(m, [t=1:NumTime], H2HP[t] * 3412.14 == (HP2PCM_H[t] + HP2H[t])/COP_H); # [BTU/hr]

            # Cooling (from kW to BTU/hr), node at heat pump cooling option
            @constraint(m, [t=1:NumTime], H2C[t] * 3412.14 == (C2PCM_C[t] + C2H[t])/COP_C); # [BTU/hr]

            # Heat pump power capacity constraint, node at heat pump heating mode
            @constraint(m, [t=1:NumTime], H2HP[t] <= Cap_HP_H); # [kW]

            # Heat pump power capacity constraint, node at heat pump cooling mode
            @constraint(m, [t=1:NumTime], H2C[t] <= Cap_HP_C); # [kW]

            # PCM heating storage balance constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime-1], InStoragePCM_H[t+1] == InStoragePCM_H[t] + δt * (HP2PCM_H[t] - PCM_H2H[t])); # [BTU]
            
            # PCM cooling storage balance constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime-1], InStoragePCM_C[t+1] == InStoragePCM_C[t] + δt * (C2PCM_C[t] - PCM_C2H[t])); # [BTU]
            
            # PCM heating storage discharging constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime], δt * PCM_H2H[t] <= InStoragePCM_H[t]); # [BTU]
            
            # PCM cooling storage discharging constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime], δt * PCM_C2H[t] <= InStoragePCM_C[t]); # [BTU]
            
            # PCM heating storage size constraint, node at PCM heating storage
            @constraint(m, [t=1:NumTime], InStoragePCM_H[t] <= PCM_H_Size); # [BTU]

            # PCM cooling storage size constraint, node at PCM cooling storage
            @constraint(m, [t=1:NumTime], InStoragePCM_C[t] <= PCM_C_Size); # [BTU]
        end
        
        # Single Operation Choice Constraints
        begin
            # If HP_OPTION = 1, heating mode on, heating only
            @constraint(m, [t=1:NumTime], H2HP[t] <= M * HP_OPTION[t]) # [kW]          

            # If HP_OPTION = 0, cooling mode on, cooling only
            @constraint(m, [t=1:NumTime], H2C[t] <= M * (1 - HP_OPTION[t])) # [kW]    
    
            # If TES_Charging_OPTION = 1, charge hot TES only
            @constraint(m, [t=1:NumTime], HP2PCM_H[t] <= M * TES_Charging_OPTION[t]) # [BTU/hr]     
    
            # If TES_Charging_OPTION = 0, charge cold TES only
            @constraint(m, [t=1:NumTime], C2PCM_C[t] <= M * (1 - TES_Charging_OPTION[t])) # [BTU/hr]

            # If TES_Discharging_OPTION = 1, discharge hot TES only
            @constraint(m, [t=1:NumTime], PCM_H2H[t] <= M * TES_Discharging_OPTION[t]) # [BTU/hr]          
    
            # If TES_Discharging_OPTION = 0, discharge cold TES only
            @constraint(m, [t=1:NumTime], PCM_C2H[t] <= M * (1 - TES_Discharging_OPTION[t])) # [BTU/hr]

            # Total TES Charging 
            @expression(m, TES_Charge[t=1:NumTime], HP2PCM_H[t] + C2PCM_C[t]); # [BTU/hr]

            # Total TES Discharging 
            @expression(m, TES_Discharge[t=1:NumTime], PCM_H2H[t] + PCM_C2H[t]); # [BTU/hr]

            # If TES_CD_OPTION = 1, charge only
            @constraint(m, [t=1:NumTime], TES_Charge[t] <= M * TES_CD_OPTION[t]) # [BTU/hr]          

            # If TES_CD_OPTION = 0, discharge only
            @constraint(m, [t=1:NumTime], TES_Discharge[t] <= M * (1 - TES_CD_OPTION[t])) # [BTU/hr] 

            # If B_OPTION = 1, charge only
            @constraint(m, [t=1:NumTime], PV2B[t] <= M * B_OPTION[t]) # [kW]          

            # If B_OPTION = 0, discharge only
            @constraint(m, [t=1:NumTime], B2H[t] <= M * (1 - B_OPTION[t])) # [kW] 
        end

    end
    ########### Solve  ##########
    optimize!(m);
    println("Finished solving big model") 
    ########### Model Results  ##########
    begin
        # Collect J_States (the states that come from Julia MPC which will be used as input for the Julia MPC next iteration)
        J_States_new = zeros(3)
        J_States_new[1] = PVGeneration[1] * PVSize # [kW] PV power generation from time step 1
        J_States_new[2] = value.(E_total[1]) # [kW] Total Plug Load from time step 1
        J_States_new[3] = value.(InStorageBattery[1]) # [kWh] The battery SOC at the beginning of time step 1

        # Collect operational commands needed as input for FMU
        Modelica_Input = 0
        if value.(TES_Discharge[1]) > 0
            Modelica_Input = 3
        elseif value.(HP2PCM_H[1]) > 0
            Modelica_Input = 1
        elseif value.(C2PCM_C[1]) > 0
            Modelica_Input = 2
        end
    end
    
    # Total penalty for loss of load
    TotalCost = objective_value(m) # [$]
    return TotalCost, J_States_new, [Modelica_Input, Modelica_Input]
end
=#