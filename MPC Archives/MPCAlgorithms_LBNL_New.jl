# MPC Algorithm
#=
This function is being called by Julia MPC (MPC_S V4.1.jl) as the actual MPC algorithm. The MPC
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

function OptimizeV4_1(iteration, steps, stepsize, Weather, M_States, J_States)
    ########## Instructions  ##########
    # current optimization model output for modelica (0 to 3)
    ENV["GRB_LICENSE_FILE"] = "C:\\Users\\Fred\\gurobi.lic"

    ########## Data Preparations  ##########  
    begin
        # Set timesteps 
        NumTime = size(Weather)[1];
        δt = steps/60;

        # Declare time series variables
        TemperatureAmbient = Weather.Temperature;
        PercentLighting = Weather.Lighting;
        PercentPlug = Weather.Plugs;
        PercentOccupied = Weather.Occupancy;
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
            Total_Load_0 = E_Modelica_0 + (E_other_0)*δt[1] # [kWh] Actual total Electrical Load at last timestep
            
            # Small Optimization Model to determine initial battery SOC using J_States
            begin
                println("Start solving small model")
                println("PVGen_0: ", PVGen_0)
                println("Total_Load_0: ", Total_Load_0)
                println("B_SOC_0: ", B_SOC_0)
                
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
                @constraint(m_0, InStorageBattery_past[2] == InStorageBattery_past[1] * (1 - δt[1] * BatteryLoss) + δt[1] * (PV2B_past * η - B2H_past)); # [kWh]
                @constraint(m_0, Total_Load_0 == PV2H_past * η_PVIV * δt[1] + B2H_past * η * δt[1]) # [kWh]

                # Battery discharging constraint, node at battery 
                @constraint(m_0, InStorageBattery_past[1] == B_SOC_0); # [kWh]

                # Battery discharging constraint, node at battery 
                @constraint(m_0, δt[1] * B2H_past <= InStorageBattery_past[1]); # [kWh]
                
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
                
                println()
                println("Battery SOC is: ", InStorageBattery_1)
                println()
                println("Hot PCM SOC is: ", InStoragePCM_H_1)
                println()
                println("Cold PCM SOC is: ", InStoragePCM_C_1)
                println()
                println("Indoor Temperature is: ", TemperatureIndoor_1, " [°F]")
                println()
                
            end
        end    
    end
    println("Start solving big model")  
    ########## Declare model  ##########
    begin
        # Define the model name and solver. In this case, model name is "m"
        begin
            m = Model(Gurobi.Optimizer)
            set_optimizer_attribute(m, "Method", 0)  # 2 corresponds to the barrier method
            set_optimizer_attribute(m, "FeasibilityTol", 1e-6)
            set_optimizer_attribute(m, "IntFeasTol", 1e-6)
            set_optimizer_attribute(m, "OptimalityTol", 1e-8)
            set_optimizer_attribute(m, "Threads", 4)
        end
    end
    ######## Decision variables ########
    begin
        @variable(m, PV2H[1:NumTime] >= 0); # [kW] electrical power transfer from PV to home (Berg)

        @variable(m, PV2G[1:NumTime] >= 0); # [kW] electrical power transfer from PV to ground (curtailment)

        @variable(m, PV2B[1:NumTime] >= 0); # [kW] electrical power transfer from PV to battery

        @variable(m, B2H[1:NumTime] >= 0); # [kW] electrical power transfer from battery to home (Berg)
        
        #=
            (HP_state[1, t] = 1) ==> Stand Alone Mode ON;
            (HP_state[2, t] = 1) ==> HP Heating Mode ON;
            (HP_state[3, t] = 1) ==> HP Cooling Mode ON;
        =#
        @variable(m, HP_state[1:3, 1:NumTime], Bin);

        @variable(m, heating[1:NumTime], Bin); # 1: Heating ON; 0: Heating OFF
        
        @variable(m, cooling[1:NumTime], Bin); # 1: Cooling ON; 0: Cooling OFF

        @variable(m, heating_start[1:NumTime], Bin); # 1: Heating should start; 0: Heating shouldn't start
        
        @variable(m, heating_end[1:NumTime], Bin); # 1: Heating should end; 0: Heating shouldn't end

        @variable(m, cooling_start[1:NumTime], Bin); # 1: Cooling should start; 0: Cooling shouldn't start
        
        @variable(m, cooling_end[1:NumTime], Bin); # 1: Cooling should end; 0: Cooling shouldn't end
        
        #=
            (PCM_state[1, t] = 1) ==> Stand Alone Mode ON;
            (PCM_state[2, t] = 1) ==> Charging PCM Hot Mode ON;
            (PCM_state[3, t] = 1) ==> Charging PCM Cold Mode ON;
            (PCM_state[4, t] = 1) ==> Discharging PCM Hot Mode ON;
            (PCM_state[5, t] = 1) ==> Discharging PCM Cold Mode ON;
        =#
        @variable(m, PCM_state[1:5, 1:NumTime], Bin); 

        @variable(m, PCM_SOC[1:2, 1:NumTime] >= -0.05); # PCM_SOC[1, t] ==> PCM Hot SOC; PCM_SOC[2, t] ==> PCM Cold SOC

        @variable(m, PCM_H_charge_rate[1:NumTime] >= 0); # [BTU/hr] thermal power transfer from HP heating to PCM hot storage

        @variable(m, PCM_C_charge_rate[1:NumTime] >= 0); # [BTU/hr] thermal power transfer from HP cooling to PCM cold storage

        @variable(m, InStorageBattery[1:NumTime] >= 0); # [kWh] Battery Remaining Charge

        @variable(m, TemperatureIndoor[1:NumTime] >= 0); # [°F] Indoor Air Temperature

        @variable(m, HP2R[t=1:NumTime] >= 0); # [BTU/hr] Heat delievered to the room from HP heating

        @variable(m, C2R[1:NumTime] >= 0); # [BTU/hr] Coolth delievered to the room from HP cooling

        @variable(m, G2H[1:NumTime] >= 0); # [kW] Loss of Load

        @variable(m. ON_Ratio[1:NumTime] >= 0)  # for 3 hour time steps, there is an option to scale the actions with time resolution of 5 mins, but now the model becomes nonlinear if they multiply.
    end
    ############ Objective Functions #############
    begin
        # Set single objective for minimizing loss of load in the current iteration of optimization horizon
        
        @objective(m, Min, sum(δt[t] * G2H[t] for t = 1:NumTime)); # [kWh]
    end
    ############# Expressions ############
    begin
        # DeltaTemp 
        @expression(m, TempDelta[t=1:NumTime], TemperatureAmbient[t] - TemperatureIndoor[t]); # [°F]

        # Power usage from lighting 
        @expression(m, E_Lighting[t=1:NumTime], PeakLighting * PercentLighting[t]); # [kW]

        # Power usage from plugs
        @expression(m, E_Plugs[t=1:NumTime], PeakPlugLoad * PercentPlug[t]); # [kW]

        # Power usage from Pump
        @expression(m, E_Pumps[t=1:NumTime], sum(PCM_state[ω, t] for ω = 2:5) * P_Pumps); # [kW]

        # Total power usage
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

        # Heat or "Coolth" thermal power delivered to the room directly from Heat Pump
        # @expression(m, HP2R[t=1:NumTime], H2HP[t] * COP_H * 3412.14 - PCM_H_Charge[t]); # [BTU/hr]
        # @expression(m, C2R[t=1:NumTime], H2C[t] * COP_C * 3412.14 - PCM_C_Charge[t]); # [BTU/hr]

        # Total Heat or "Coolth" thermal power delivered to the room (directly from Heat Pump or PCM discharge)
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
        (δt[t]/TC)*(Q_Others[t] + Q_Envelope[t] + Q_Rad[t] - Q_RadCool[t] + Heat_delivered[t] - Cool_delivered[t])); # [°F]
        
        # PV power balance constraint, node at PV
        @constraint(m, [t=1:NumTime], PVGeneration[t] * PVSize == PV2B[t] + PV2H[t] + PV2G[t]); # [kW]
    
        # House power load balance constraint, node at house, battery efficiency modeled
        @constraint(m, [t=1:NumTime], E_total[t] + H2HP[t] + H2C[t] == PV2H[t] * η_PVIV + B2H[t] * η + G2H[t]); # [kW]

        # Battery storage balance constraint, node at battery, battery leakage modeled, battery efficiency modeled
        @constraint(m, [t=1:NumTime-1], InStorageBattery[t+1] == InStorageBattery[t] * (1 - δt[t] * BatteryLoss) +  δt[t] * (PV2B[t] * η - B2H[t])); # [kWh]

        # Battery discharging constraint, node at battery 
        @constraint(m, [t=1:NumTime-1], δt[t] * B2H[t] <= InStorageBattery[t]); # [kWh]
        
        # Battery power inverter constraint, node at battery (inverter power constraint) 
        @constraint(m, [t=1:NumTime], B2H[t] + PV2B[t] <= InverterSize); # [kW]
        
        # Battery storage size constraint, node at battery
        @constraint(m, [t=1:NumTime], InStorageBattery[t] <= BatterySize); # [kWh]
        
        # Battery storage max discharge constraint, node at battery, always at least 20% full 
        @constraint(m, [t=1:NumTime], InStorageBattery[t] >= BatterySize * (1-MaxDischarge)); # [kWh]
       
        # Heating and Cooling band temperature controls
        begin 
            # Basically, all the following constraints in this section are made to achieve the following hysteresis logic modeled in Modelica:
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
                # Given current indoor temperature, decide whether heating should start or end
                # If indoor temperature <= (zone_temp_heating_setpoint - zone_temp_setpoint_delta), heating should start and not end;
                # If indoor temperature >= (zone_temp_heating_setpoint + zone_temp_setpoint_delta), heating should end and not start;
                # If indoor temperature is in between, then heating shouldn't start nor end, but keeping the same mode from last timestep;
                @constraint(m, [t=1:NumTime], TemperatureIndoor[t] - (zone_temp_heating_setpoint - zone_temp_setpoint_delta) <= (1 - heating_start[t]) * M)
                @constraint(m, [t=1:NumTime], (zone_temp_heating_setpoint - zone_temp_setpoint_delta) - TemperatureIndoor[t] <= heating_start[t] * M)
                @constraint(m, [t=1:NumTime], TemperatureIndoor[t] - (zone_temp_heating_setpoint + zone_temp_setpoint_delta) >= (1 - heating_end[t]) * -M)
                @constraint(m, [t=1:NumTime], (zone_temp_heating_setpoint + zone_temp_setpoint_delta) - TemperatureIndoor[t] >= heating_end[t] * -M) 

                # Given heating_start, heating_end, and the last timestep's Heating Mode, decide whether Heating Mode should be on
                # Ensure that Heating Mode is ON (heating[t] = 1) if heating_start[t] == 1 and heating_end[t] == 0;
                # Ensure that Heating Mode is OFF (heating[t] = 0) if heating_start[t] == 0 and heating_end[t] == 1;
                # Ensure that Heating Mode is ON (heating[t] = 1) if heating_start[t] == 0 and heating_end[t] == 0, but heating[t-1] = 1;
                # Ensure that Heating Mode is OFF (heating[t] = 0) if heating_start[t] == 0 and heating_end[t] == 0, but heating[t-1] = 0;            
                begin
                    @constraint(m, [t=1:NumTime], heating[t] >= heating_start[t] - heating_end[t]); # [Bin]
                                
                    # Additional constraints to ensure the logic is maintained over time
                    @constraint(m, heating[1] <= heating_0 + heating_start[1]); # [Bin]
                    @constraint(m, heating[1] <= 1 - heating_end[1]); # [Bin]
                    @constraint(m, heating[1] >= heating_0 - heating_end[1]); # [Bin]

                    @constraint(m, [t=2:NumTime], heating[t] <= heating[t-1] + heating_start[t]); # [Bin]
                    @constraint(m, [t=2:NumTime], heating[t] <= 1 - heating_end[t]); # [Bin]
                    @constraint(m, [t=2:NumTime], heating[t] >= heating[t-1] - heating_end[t]); # [Bin]
                end
            end
            
            # Cooling Logic
            begin
                # Given current indoor temperature, decide whether cooling should start or end
                # If indoor temperature >= (zone_temp_cooling_setpoint + zone_temp_setpoint_delta), cooling should start and not end;
                # If indoor temperature <= (zone_temp_cooling_setpoint - zone_temp_setpoint_delta), cooling should end and not start;
                # If indoor temperature is in between, then cooling shouldn't start nor end, but keeping the same mode from last timestep;
                @constraint(m, [t=1:NumTime], TemperatureIndoor[t] - (zone_temp_cooling_setpoint + zone_temp_setpoint_delta) <= (1 - cooling_start[t]) * M)
                @constraint(m, [t=1:NumTime], (zone_temp_cooling_setpoint + zone_temp_setpoint_delta) - TemperatureIndoor[t] <= cooling_start[t] * M)
                @constraint(m, [t=1:NumTime], TemperatureIndoor[t] - (zone_temp_cooling_setpoint - zone_temp_setpoint_delta) >= (1 - cooling_end[t]) * -M)
                @constraint(m, [t=1:NumTime], (zone_temp_cooling_setpoint - zone_temp_setpoint_delta) - TemperatureIndoor[t] >= cooling_end[t] * -M)

                # Given cooling_start, cooling_end, and the last timestep's Cooling Mode, decide whether Cooling Mode should be on
                # Ensure that Cooling Mode is ON (cooling[t] = 1) if cooling_start[t] == 1 and cooling_end[t] == 0;
                # Ensure that Cooling Mode is OFF (cooling[t] = 0) if cooling_start[t] == 0 and cooling_end[t] == 1;
                # Ensure that Cooling Mode is ON (cooling[t] = 1) if cooling_start[t] == 0 and cooling_end[t] == 0, but cooling_end[t-1] = 1;
                # Ensure that Cooling Mode is OFF (cooling[t] = 0) if cooling_start[t] == 0 and cooling_end[t] == 0, but cooling_end[t-1] = 0;            
                
                begin
                    @constraint(m, [t=1:NumTime], cooling[t] >= cooling_start[t] - cooling_end[t]); # [Bin]
                
                    # Additional constraints to ensure the logic is maintained over time
                    @constraint(m, cooling[1] <= cooling_0 + cooling_start[1]); # [Bin]
                    @constraint(m, cooling[1] <= 1 - cooling_end[1]); # [Bin]
                    @constraint(m, cooling[1] >= cooling_0 - cooling_end[1]); # [Bin]
    
                    @constraint(m, [t=2:NumTime], cooling[t] <= cooling[t-1] + cooling_start[t]); # [Bin]
                    @constraint(m, [t=2:NumTime], cooling[t] <= 1 - cooling_end[t]); # [Bin]
                    @constraint(m, [t=2:NumTime], cooling[t] >= cooling[t-1] - cooling_end[t]); # [Bin]
                end
            end
            
            # If Heating Mode is ON, ensure that at least one between HP and PCM H Discharge is heating the room
            @constraint(m, [t=1:NumTime], HP_state[2, t] + PCM_state[4, t] <= 2 * heating[t]); # [Bin]            
            @constraint(m, [t=1:NumTime], HP_state[2, t] + PCM_state[4, t] >= heating[t]); # [Bin]
            
            # If Cooling Mode is ON, ensure that at least one between HP and PCM C Discharge is cooling the room
            @constraint(m, [t=1:NumTime], HP_state[3, t] + PCM_state[5, t] <= 2 * cooling[t]); # [Bin]
            @constraint(m, [t=1:NumTime], HP_state[3, t] + PCM_state[5, t] >= cooling[t]); # [Bin]
            
            # Heat Pump thermal power balance constraints, node at Heat Pump
            # @constraint(m, [t=1:NumTime], HP2R[t] >= 0); # [BTU/hr]
            # @constraint(m, [t=1:NumTime], C2R[t] >= 0); # [BTU/hr]

            @constraint(m, [t=1:NumTime], HP2R[t] == H2HP[t] * COP_H * 3412.14 - PCM_H_Charge[t]); # [BTU/hr]
            @constraint(m, [t=1:NumTime], C2R[t] == H2C[t] * COP_C * 3412.14 - PCM_C_Charge[t]); # [BTU/hr]
        end

        # PCM heating storage initialization constraint, node at PCM heating storage
        @constraint(m, PCM_SOC[1, 1] == InStoragePCM_H_1); # [1]

        # PCM cooling storage initialization constraint, node at PCM cooling storage
        @constraint(m, PCM_SOC[2, 1] == InStoragePCM_C_1); # [1]

        # PCM heating storage balance constraint, node at PCM heating storage
        @constraint(m, [t=1:NumTime-1], PCM_SOC[1, t+1] == PCM_SOC[1, t] + (δt[t] * (PCM_H_Charge[t] - PCM_H_Discharge[t]))/PCM_H_Size); # [1]
        
        # PCM cooling storage balance constraint, node at PCM cooling storage
        @constraint(m, [t=1:NumTime-1], PCM_SOC[2, t+1] == PCM_SOC[2, t] + (δt[t] * (PCM_C_Charge[t] - PCM_C_Discharge[t]))/PCM_C_Size); # [1]
        
        # PCM heating and cooling storage state of charge constraint, node at PCM storage
        @constraint(m, [ω = 1:2, t=2:NumTime], PCM_SOC[ω, t] <= 1.05); # [1]

        # @constraint(m, [t=1:NumTime], PCM_H_charge_rate[t] <= 3412.14*0.5); # [BTU/hr] thermal power transfer from HP heating to PCM hot storage

        # @constraint(m, [t=1:NumTime], PCM_C_charge_rate[t] <= 3412.14*0.5); # [BTU/hr] thermal power transfer from HP cooling to PCM cold storage
    end

    ########### Solve  ##########
    optimize!(m);
    # println(m)
    # Check the solver status
    status = JuMP.termination_status(m)
    println(status)
    println("Finished solving big model") 
    
    Didnt_Solve = 1;
    ########### Model Results  ##########
    begin
        if status == MOI.OPTIMAL
            # Collect J_States (Julia States: the states that come from Julia MPC which will be used as input for the Julia MPC at next iteration)
            J_States_new = zeros(3)

            J_States_new[1] = PVGeneration[1] * PVSize # [kW] PV power generation from current timestep
            J_States_new[2] = value.(E_total[1]) # [kW] Total Load from time step 1
            J_States_new[3] = value.(InStorageBattery[1]) # [kWh] The battery SOC at the beginning of time step 1

            println(J_States_new)

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

            # Total penalty for loss of load
            CurrentCost = δt[1] * value.(G2H[1]) # [kWh]

            Didnt_Solve = 0;
        else
            # If the model is not feasible, return default values
            println("Model is infeasible or solver encountered an issue.")
            J_States_new = [PVGeneration[1] * PVSize, (PeakLighting * PercentLighting[1] + PeakPlugLoad * PercentPlug[1] + P_Controls + P_AHU), InStorageBattery_1] # Default values for J_States
            Modelica_Input = 0              # Default value for Modelica_Input
            CurrentCost = 0.0               # Default cost
        end
    end

    # Return results
    return CurrentCost, J_States_new, [Modelica_Input, Modelica_Input], Didnt_Solve
end
