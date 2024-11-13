using JuMP
using Gurobi
ENV["GRB_LICENSE_FILE"] = "C:\\Users\\Fred\\.julia\\environments\\v1.7\\gurobi.lic" # Fred's license
# Create a model
model = Model(Gurobi.Optimizer)

# Define variables
@variable(model, x >= 0)
@variable(model, y >= 0)

# Define some conditions based on external data or logic
external_data1 = 8  # Example external data for condition 1
external_data2 = 9  # Example external data for condition 2
a = 9

# Define conditions
function condition1(data)
    return data < 10
end

function condition2(data)
    return data > 5
end

# Add rule-based constraints
function add_rule_based_constraints(model, data1, data2, a)
    if condition1(data1)
        a = a + 1
        @constraint(model, x <= a)
    end
    if condition2(data2)
        @constraint(model, y >= 5)
    end
end

# Add constraints based on conditions
add_rule_based_constraints(model, external_data1, external_data2, a)

# Objective function
@objective(model, Min, x + y)

# Solve the model
optimize!(model)

# Get the solution
println("x = ", value(x))
println("y = ", value(y))



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
@variable(model, switching[1:T-1], Bin)  # Binary variable for state switching
@variable(model, thermal_storage_state[1:T, 1:4], Bin)  # Four states: do nothing, charge hot, charge cold, discharge

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


# Heating start and end logic
@constraint(model, heating_start[t] >= (zone_temp_heating_setpoint - zone_temp_setpoint_delta - indoor_temp[t]) / M)
@constraint(model, heating_start[t] <= 1 - heating_end[t-1])
@constraint(model, heating_end[t] >= (indoor_temp[t] - (zone_temp_heating_setpoint + zone_temp_setpoint_delta)) / M)
@constraint(model, heating_end[t] <= heating[t-1])





# Define the constraints based on the rules provided
for t in T
    @constraint(model, heating_start[t] == (Ti[t] <= zone_temp_heating_setpoint - zone_temp_setpoint_delta ? 1 : 0))
    @constraint(model, heating_end[t] == (Ti[t] >= zone_temp_heating_setpoint + zone_temp_setpoint_delta ? 1 : 0))
    @constraint(model, heating[t] <= heating_start[t])
    @constraint(model, heating[t] <= 1 - heating_end[t])
end

# Define the constraints based on the rules provided
for t in T
    @constraint(model, cooling_start[t] == (Ti[t] >= zone_temp_cooling_setpoint + zone_temp_setpoint_delta ? 1 : 0))
    @constraint(model, cooling_end[t] == (Ti[t] <= zone_temp_cooling_setpoint - zone_temp_setpoint_delta ? 1 : 0))
    @constraint(model, cooling[t] <= cooling_start[t])
    @constraint(model, cooling[t] <= 1 - cooling_end[t])
end

