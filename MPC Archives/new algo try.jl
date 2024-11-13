using JuMP, GLPK  # Use GLPK or Gurobi if you have a license

# Model
model = Model(GLPK.Optimizer)

# Parameters
T = 168  # Number of hours in a week
zone_temp_cooling_setpoint = 273.15 + 22
zone_temp_heating_setpoint = 273.15 + 20
zone_temp_setpoint_delta = 0.5
max_switches = 10  # Maximum allowed switches in the time span


# Variables
@variable(model, pv_gen[1:T] >= 0)
@variable(model, battery_charge[1:T] >= 0)
@variable(model, battery_discharge[1:T] >= 0)
@variable(model, indoor_temp[1:T] >= zone_temp_heating_setpoint - zone_temp_setpoint_delta, upper_bound=zone_temp_cooling_setpoint + zone_temp_setpoint_delta)
@variable(model, heating[1:T], Bin)
@variable(model, cooling[1:T], Bin)
@variable(model, deadband[1:T], Bin)
@variable(model, thermal_storage_state[1:T, 1:4], Bin)  # Four states: do nothing, charge hot, charge cold, discharge
@variable(model, switching[1:T-1], Bin)  # Binary variable for state switching

# Initial conditions
@constraint(model, indoor_temp[1] == 273.15 + 20)  # Assuming an initial indoor temperature
@constraint(model, heating[1] + cooling[1] + deadband[1] == 1)  # Initial mode constraint

# Heat pump mode constraints
for t in 2:T
    @constraint(model, heating[t] + cooling[t] + deadband[t] == 1)
    
    # Heating mode rules
    @constraint(model, indoor_temp[t] <= zone_temp_heating_setpoint + zone_temp_setpoint_delta + 1000 * (1 - heating[t]))
    @constraint(model, indoor_temp[t] >= zone_temp_heating_setpoint - zone_temp_setpoint_delta - 1000 * (1 - heating[t]))
    
    # Cooling mode rules
    @constraint(model, indoor_temp[t] >= zone_temp_cooling_setpoint - zone_temp_setpoint_delta - 1000 * (1 - cooling[t]))
    @constraint(model, indoor_temp[t] <= zone_temp_cooling_setpoint + zone_temp_setpoint_delta + 1000 * (1 - cooling[t]))
    
    # Deadband mode rules
    @constraint(model, indoor_temp[t] <= zone_temp_heating_setpoint + zone_temp_setpoint_delta + 1000 * (1 - deadband[t]))
    @constraint(model, indoor_temp[t] >= zone_temp_cooling_setpoint - zone_temp_setpoint_delta - 1000 * (1 - deadband[t]))
    
    # Switching logic
    @constraint(model, switching[t-1] >= heating[t] - heating[t-1])
    @constraint(model, switching[t-1] >= heating[t-1] - heating[t])
    @constraint(model, switching[t-1] >= cooling[t] - cooling[t-1])
    @constraint(model, switching[t-1] >= cooling[t-1] - cooling[t])
    @constraint(model, switching[t-1] >= deadband[t] - deadband[t-1])
    @constraint(model, switching[t-1] >= deadband[t-1] - deadband[t])
end

# Limit the number of switches
@constraint(model, sum(switching[t] for t in 1:T-1) <= max_switches)

# Simplified energy balance constraints (you will need to refine these based on your actual system)
for t in 1:T
    # Example: balance for battery and PV generation
    @constraint(model, battery_charge[t] + battery_discharge[t] == pv_gen[t] - thermal_storage_state[t, 2] * 10 + thermal_storage_state[t, 4] * 10)
    
    # Add more detailed energy balances for heat pump, thermal storage, etc.
end

# Objective function: Minimize total energy cost (simplified)
@objective(model, Min, sum(battery_discharge[t] for t in 1:T))

# Solve the model
optimize!(model)

# Print results
println("PV Generation: ", value.(pv_gen))
println("Battery Charge: ", value.(battery_charge))
println("Battery Discharge: ", value.(battery_discharge))
println("Indoor Temperature: ", value.(indoor_temp))
println("Heating Mode: ", value.(heating))
println("Cooling Mode: ", value.(cooling))
println("Deadband Mode: ", value.(deadband))
println("Switching: ", value.(switching))
println("Thermal Storage State: ", value.(thermal_storage_state))


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

# Old Method with Ifs
begin            
    # Heating Logic
    for t = 1:NumTime
        if TemperatureIndoor[t] <= zone_temp_heating_setpoint - zone_temp_setpoint_delta
            @constraint(m, heating_start[t] == 1)
            @constraint(m, heating_end[t] == 0)
        elseif TemperatureIndoor[t] >= zone_temp_heating_setpoint + zone_temp_setpoint_delta
            @constraint(m, heating_start[t] == 0)
            @constraint(m, heating_end[t] == 1)
        else
            @constraint(m, heating_start[t] == 0)
            @constraint(m, heating_end[t] == 0)
        end

        # Ensure that heating[t] is 1 if heating_start[t] == 1 and heating_end[t] == 0
        @constraint(m, heating[t] >= heating_start[t] - heating_end[t])
        
        # Additional constraints to ensure the logic is maintained over time
        if t == 1
            @constraint(m, heating[t] <= heating_start[t])
            @constraint(m, heating[t] <= 1 - heating_end[t])
        else
            @constraint(m, heating[t] <= heating[t-1] + heating_start[t])
            @constraint(m, heating[t] <= 1 - heating_end[t])
            @constraint(m, heating[t] >= heating[t-1] - heating_end[t])
        end
    end
    
    # Cooling Logic
    for t = 1:NumTime
        if TemperatureIndoor[t] >= zone_temp_heating_setpoint + zone_temp_setpoint_delta
            @constraint(m, cooling_start[t] == 1)
            @constraint(m, cooling_end[t] == 0)
        elseif TemperatureIndoor[t] <= zone_temp_heating_setpoint - zone_temp_setpoint_delta
            @constraint(m, cooling_start[t] == 0)
            @constraint(m, cooling_end[t] == 1)
        else
            @constraint(m, cooling_start[t] == 0)
            @constraint(m, cooling_end[t] == 0)
        end
        
        # Ensure that cooling[t] is 1 if cooling_start[t] == 1 and cooling_end[t] == 0
        @constraint(m, cooling[t] >= cooling_start[t] - cooling_end[t])
        
        # Additional constraints to ensure the logic is maintained over time
        if t == 1
            @constraint(m, cooling[t] <= cooling_start[t])
            @constraint(m, cooling[t] <= 1 - cooling_end[t])
        else
            @constraint(m, cooling[t] <= cooling[t-1] + cooling_start[t])
            @constraint(m, cooling[t] <= 1 - cooling_end[t])
            @constraint(m, cooling[t] >= cooling[t-1] - cooling_end[t])
        end
    end
end

    
    


PCM_State[2, t] == 0
PCM_State[3, t] == 0


# PCM SOC Logic
begin            
    # Heating Logic
    for t = 1:NumTime
        if SOC[t] <= 0.50
            @constraint(m, not_full[t] == 1)
            @constraint(m, full[t] == 0)
        elseif SOC[t] >= 0.98
            @constraint(m, not_full[t] == 0)
            @constraint(m, full[t] == 1)
        else
            @constraint(m, not_full[t] == 0)
            @constraint(m, full[t] == 0)
        end

        # Ensure that heating[t] is 1 if heating_start[t] == 1 and heating_end[t] == 0
        @constraint(m, chargable[t] >= not_full[t] - full[t])
        
        # Additional constraints to ensure the logic is maintained over time
        if t == 1
            @constraint(m, chargable[t] <= not_full[t])
            @constraint(m, chargable[t] <= 1 - full[t])
        else
            @constraint(m, chargable[t] <= chargable[t-1] + not_full[t])
            @constraint(m, chargable[t] <= 1 - full[t])
            @constraint(m, chargable[t] >= chargable[t-1] - full[t])
        end
    end    
end

# Big M
for t in T
    @constraint(model, heating_start[t] <= (TemperatureIndoor[t] <= zone_temp_heating_setpoint - zone_temp_setpoint_delta ? 1 : 0) * M)
    @constraint(model, heating_end[t] <= (TemperatureIndoor[t] >= zone_temp_heating_setpoint + zone_temp_setpoint_delta ? 1 : 0) * M)
    @constraint(model, heating_start[t] >= (TemperatureIndoor[t] <= zone_temp_heating_setpoint - zone_temp_setpoint_delta ? 1 : 0))
    @constraint(model, heating_end[t] >= (TemperatureIndoor[t] >= zone_temp_heating_setpoint + zone_temp_setpoint_delta ? 1 : 0))

    @constraint(model, cooling_start[t] <= (TemperatureIndoor[t] >= zone_temp_heating_setpoint + zone_temp_setpoint_delta ? 1 : 0) * M)
    @constraint(model, cooling_end[t] <= (TemperatureIndoor[t] <= zone_temp_heating_setpoint - zone_temp_setpoint_delta ? 1 : 0) * M)
    @constraint(model, cooling_start[t] >= (TemperatureIndoor[t] >= zone_temp_heating_setpoint + zone_temp_setpoint_delta ? 1 : 0))
    @constraint(model, cooling_end[t] >= (TemperatureIndoor[t] <= zone_temp_heating_setpoint - zone_temp_setpoint_delta ? 1 : 0))
end

# PCM SOC Logic
begin            
    for t = 1:NumTime
        if SOC[t] <= 0.50
            @constraint(m, not_full[t] == 1)
            @constraint(m, full[t] == 0)
        elseif SOC[t] >= 0.98
            @constraint(m, not_full[t] == 0)
            @constraint(m, full[t] == 1)
        else
            @constraint(m, not_full[t] == 0)
            @constraint(m, full[t] == 0)
        end

        # Ensure that heating[t] is 1 if heating_start[t] == 1 and heating_end[t] == 0
        @constraint(m, chargable[t] >= not_full[t] - full[t])
        
        # Additional constraints to ensure the logic is maintained over time
        if t == 1
            @constraint(m, chargable[t] <= not_full[t])
            @constraint(m, chargable[t] <= 1 - full[t])
        else
            @constraint(m, chargable[t] <= chargable[t-1] + not_full[t])
            @constraint(m, chargable[t] <= 1 - full[t])
            @constraint(m, chargable[t] >= chargable[t-1] - full[t])
        end
    end    
end


# VERY beginning
# Additional constraints to ensure the logic is maintained over time
@constraint(m, heating[1] <= heating_start[1]); # this is at the very very beginning
@constraint(m, heating[1] <= 1 - heating_end[1]);
# Additional constraints to ensure the logic is maintained over time
@constraint(m, cooling[1] <= cooling_start[1])
@constraint(m, cooling[1] <= 1 - cooling_end[1])
# same for pcm
cooling_0, heating_0, chargable_0[ω] # get these


# these can be constant
HP_power_H # is this the standard HP power rate only for home or for both charging PCM H and home? how is this related to PCM_H_charge_rate?
HP_power_C
PCM_H_discharge_rate
PCM_C_discharge_rate

# variable COP H and C
# Heating and Cooling Constraints (Needs Remodeling)
begin
    # DeltaTemp 
    # @expression(m, TempDelta_C[t=1:NumTime], (TemperatureAmbient[t] * 1.8 + 32) - (TemperatureIndoor[t] * 1.8 + 32)); # [°C]

    # Detailed COP of heating and cooling (requires nonlinear optimization solver, although we can hold indoor temperature as a constant to keep the model linear)
    # COP (Heating Home)
    # @expression(m, COP_Heating_H[t=1:NumTime], HP_a + HP_b * (-TempDelta_C[t]) + HP_c * TempDelta[t]^2);
    
    # COP (Heating PCM H)
    # @expression(m, COP_Heating_PCM[t=1:NumTime], HP_a + HP_b * (48 - Weather[t, 9]) + HP_c * (48 - Weather[t, 9])^2); 

    # COP (Cooling Home)
    # @expression(m, COP_Cooling_H[t=1:NumTime], HP_a + HP_b * (TempDelta_C[t]) + HP_c * TempDelta[t]^2);
    
    # COP (Cooling PCM C)
    # @expression(m, COP_Cooling_PCM[t=1:NumTime], HP_a + HP_b * (Weather[t, 9] - 11) + HP_c * (Weather[t, 9] - 11)^2); 

    #=
    # Detailed COP of heating and cooling (requires nonlinear optimization solver)
    # Heating (from kWh to BTU), node at heat pump heating option
    @constraint(m, [t=1:NumTime], H2HP[t] * 3412.14 == HP2PCM_H[t]/COP_Heating_PCM[t] + HP2H[t]/COP_Heating_H[t]);
    
    # Cooling (from kWh to BTU), node at heat pump cooling option
    @constraint(m, [t=1:NumTime], H2C[t] * 3412.14 == C2PCM_C[t]/COP_Cooling_PCM[t] + C2H[t]/COP_Cooling_H[t]);
    =#
end

# can heat pump simutaneously heat ambient and charge PCM H, whats the max allowed rate?
# I didn't include the capacity of HP for this reason
# Heat pump power capacity constraint, node at heat pump heating mode
@constraint(m, [t=1:NumTime], H2HP[t] <= Cap_HP_H); # [kW]

# Heat pump power capacity constraint, node at heat pump cooling mode
@constraint(m, [t=1:NumTime], H2C[t] <= Cap_HP_C); # [kW]

@variable(m, H2HP[1:NumTime] >= 0); # [kW] electrical power transfer from home (Berg) to heat pump heating mode

@variable(m, HP2H[1:NumTime] >= 0); # [BTU/hr] heating power transfer from heat pump heating mode to home (Berg)

@variable(m, H2C[1:NumTime] >= 0); # [kW] electrical power transfer from home (Berg) to heat pump cooling unit

@variable(m, C2H[1:NumTime] >= 0); # [BTU/hr] cooling power transfer from heat pump cooling mode to home (Berg)

@variable(m, InStoragePCM_H[1:NumTime] >= 0); # [BTU] PCM Heating Remaining Charge

@variable(m, InStoragePCM_C[1:NumTime] >= 0); # [BTU] PCM Cooling Remaining Charge

# Need to model pump and fan electricity consumption, need constants from LBNL

# Need to model PCM SOC curve, need constants from LBNL

# PCM "fully charged" state from Modelica