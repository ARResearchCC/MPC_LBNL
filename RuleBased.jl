using Pkg

using Gurobi
using JuMP
using DataFrames

include("Input_Parameters.jl")
include("User_Input.jl")

# input_df needs ambient temperature as a column DONE
# Divide QPassive into heating_load and cooling_load DONE
# if knowing the precise amount of heating/cooling power [kW] of each individual heating/cooling call, that would be great.
# input previous decision 0-3, i = 1, or i > 1
# PV input is already in kWh, and make sure battery SOC gets converted to kWh

function RuleBased(iteration, input_df, M_States, J_States, previous_decision)

    # Set timesteps 
    NumTime = size(input_df)[1];

    # Length of the M_States dataframe
    msize = size(M_States)[1];

    PVGeneration = input_df[1, :PV] # [kW]
    E_load = input_df[1, :E_load] # [kW]

    PCM_C_temp = M_States[end, :"PCM_Cold_Temp"]
    PCM_H_temp = M_States[end, :"PCM_Hot_Temp"]

    PCM_C_Energy, PCM_H_Energy = PCM_available_energy(PCM_C_temp, PCM_H_temp) # also set a copy as initial variable

    # T_HP = 20 
    # compressor_speed = 0
    Ta = input_df[1, :Ta]
    option = "NULL"

    # These two variables mean that whether the system is in heating mode or cooling mode at the start of the current iteration.
    # heating_0 = 0
    # cooling_0 = 0

    time_intervals = zeros(msize-1); # [hr]
    
    heating_0 = 0
    cooling_0 = 0

    if M_States[msize, :"FCU Mode"] == 1
        heating_0 = 1
    end
    if M_States[msize, :"FCU Mode"] == -1
        cooling_0 = 1
    end

    for i = 1:msize-1
        time_intervals[i] = (M_States[i+1, :"Time"] - M_States[i, :"Time"])/3600 # [hr]
    end

    E_Pump_1_0 = sum(((M_States[i+1, :"Pump 1 Electric Power"] + M_States[i, :"Pump 1 Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]

    E_Pump_2_0 = sum(((M_States[i+1, :"Pump 2 Electric Power"] + M_States[i, :"Pump 2 Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
    
    E_HP_0 = sum(((M_States[i+1, :"HP Electric Power"] + M_States[i, :"HP Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
    
    E_Fan_0 = sum(((M_States[i+1, :"Fan Coil Fan Electric Power"] + M_States[i, :"Fan Coil Fan Electric Power"])/2)*time_intervals[i] for i = 1:(msize-1))/1000 # [kWh]
    
    # Julia States
    PVGen_0 = J_States[1, 1] # [kW] PV Generation during last timestep (Actual)
    E_other_0 = J_States[1, 2] # [kW] Other Electrical Load (from schedule) during last timestep (Actual)
    B_SOC_0 = J_States[1, 3] # [kWh] Battery SOC at the start of last timestep (Actual)
    
    E_Modelica_0 = (E_Pump_1_0 + E_Pump_2_0 + E_HP_0 + E_Fan_0)/δt # [kW] Actual Electrical Load from HP system at last timestep
    Total_Load_0 = E_Modelica_0 + E_other_0 # [kW] Actual total Electrical Load during last timestep

    T_HP = M_States[msize, :"Heat Pump Return Water Temperature"]
    compressor_speed = M_States[msize, :"Heat Pump Compressor Speed"]

    # Using the initial battery SOC at the start of the last timestep, PV generation during last timestep, the actual total load during last timestep,
    # we are able to solve for the curtailment OR loss of load happend during the last timestep, and the end battery SOC at the last timestep, which is the initial battery SOC at the start of the current timestep.
    SOC, curtailment, loss_of_load = Get_battery_SOC(PVGen_0, Total_Load_0, B_SOC_0) 

    SOC_initial = copy(SOC)
    PCM_C_Energy_initial = copy(PCM_C_Energy)
    PCM_H_Energy_initial = copy(PCM_H_Energy)
    
    Decision = 999
    Scenario = 999
    loss_of_load_E = 0 # [kW] unused
    loss_of_load_Q = 0 # [kW] unused
    # curtailment = 0 # [kW] unused

    PV2B = 0
    PV2H = 0
    PV2G = 0
    B2H = 0
    G2H = 0
    Available_PV = 0 # [kW]

    PCM_H_charge = 0
    PCM_C_charge = 0

    PCM_H_discharge = 0
    PCM_C_discharge = 0


    # thermal power: we assume that a standard heating or cooling call requires 3 kW of thermal power
    if heating_0 == 1  
        heating_power = 3 # [kW]
    elseif cooling_0 == 1
        cooling_power = 3 # [kW]
    end

    if PVGeneration == 0 # Night Time
        SOC, loss_of_load = discharge_battery(SOC, E_load) # Discharge battery to meet electrical load first
        if loss_of_load > 0 # Loss of load scenario
            loss_of_load_E = loss_of_load_E + loss_of_load
            if heating_0 == 1
                option = "H"
                PCM_H_Energy, owed_power = discharge_PCM(option, PCM_H_Energy, heating_power)
                if owed_power > 0
                    loss_of_load_Q = owed_power
                    Decision = 3
                    Scenario = 4
                else
                    Decision = 3
                    Scenario = 2
                end
            elseif cooling_0 == 1
                option = "C"
                PCM_C_Energy, owed_power = discharge_PCM(option, PCM_C_Energy, cooling_power)
                if owed_power > 0
                    loss_of_load_Q = owed_power
                    Decision = 3
                    Scenario = 8
                else
                    Decision = 3
                    Scenario = 6
                end
            else
                Decision = 0
                Scenario = 0
            end
        else
            if heating_0 == 1
                option = "H"
                favorable = check_favorable_COP(option, input_df, PCM_H_Energy)
                if favorable === true # heating call ON, HP heating favorable, electricity to thermal conversion favorable, prioritize using electricity
                    SOC, extra_thermal_power, owed_thermal_power = discharge_battery_to_HP(option, option, SOC, max_heating_power, 0, Ta, T_HP, compressor_speed) # discharge battery to use HP
                    if owed_thermal_power > 0 # heating call not met by HP powered by battery discharge
                        PCM_H_Energy, owed_thermal_power = discharge_PCM(option, PCM_H_Energy, owed_thermal_power) # discharge PCM H
                        if owed_thermal_power > 0
                            loss_of_load_Q = owed_thermal_power
                        end
                        Decision = 3
                        Scenario = 4
                    else # heating call met by HP powered by battery discharge, extra heat produced, charge PCM H (hybrid charge mode)
                        PCM_H_Energy, left_over_thermal = charge_PCM_simple(option, PCM_H_Energy, extra_thermal_power)
                        Decision = 1
                        Scenario = 3
                    end
                else # heating call ON, HP heating unfavorable (or not urgent enough, therefore unneeded), prioritize using thermal
                    PCM_H_Energy, owed_thermal_power = discharge_PCM(option, PCM_H_Energy, heating_power) # discharge PCM H
                    if owed_thermal_power > 0 # heating call not met by PCM H discharge alone
                        SOC, extra_thermal_power, owed_thermal_power = discharge_battery_to_HP(option, option, SOC, owed_thermal_power, 0, Ta, T_HP, compressor_speed) # this is not great, as HP should be operating in standard power
                        if owed_thermal_power > 0
                            loss_of_load_Q = owed_thermal_power
                        end
                        Decision = 3
                        Scenario = 4
                    else # heating call met by discharging PCM H alone
                        Decision = 3
                        Scenario = 2
                    end
                end    
            elseif cooling_0 == 1 # need optimization
                option = "C"
                favorable = check_favorable_COP(option, input_df, PCM_C_Energy)
                if favorable === true # cooling call ON, HP cooling favorable, electricity to thermal conversion favorable, prioritize using electricity
                    SOC, extra_thermal_power, owed_thermal_power = discharge_battery_to_HP(option, option, SOC, max_cooling_power, 0, Ta, T_HP, compressor_speed) # discharge battery to use HP
                    if owed_thermal_power > 0 # cooling call not met by HP powered by battery discharge
                        PCM_C_Energy, owed_thermal_power = discharge_PCM(option, PCM_C_Energy, owed_thermal_power) # discharge PCM C
                        if owed_thermal_power > 0
                            loss_of_load_Q = owed_thermal_power
                        end
                        Decision = 3
                        Scenario = 8
                    else # cooling call met by HP powered by battery discharge, extra coolth produced, charge PCM C (hybrid charge mode)
                        PCM_C_Energy, left_over_thermal = charge_PCM_simple(option, PCM_C_Energy, extra_thermal_power)
                        Decision = 2
                        Scenario = 7
                    end
                else # cooling call ON, HP cooling unfavorable (or not urgent enough, therefore unneeded), prioritize using thermal
                    PCM_C_Energy, owed_thermal_power = discharge_PCM(option, PCM_C_Energy, heating_power) # discharge PCM C
                    if owed_thermal_power > 0 # cooling call not met by PCM C discharge alone
                        SOC, extra_thermal_power, owed_thermal_power = discharge_battery_to_HP(option, option, SOC, owed_thermal_power, 0, Ta, T_HP, compressor_speed) # this is not great, as HP should be operating in standard power
                        if owed_thermal_power > 0
                            loss_of_load_Q = owed_thermal_power
                        end
                        Decision = 3
                        Scenario = 8
                    else # cooling call met by discharging PCM C alone
                        Decision = 3
                        Scenario = 6
                    end 
                end
            else # no heating or cooling Call
                urgency_H, urgency_C = Get_urgency_score(input_df, M_States)
                if urgency_H >= 0.7 # heating is urgent
                    option = "H"
                    favorable = check_favorable_COP(option, input_df, PCM_H_Energy)
                    if favorable === true # charge PCM H
                        SOC, extra_thermal_power, owed_thermal_power = discharge_battery_to_HP('N', option, SOC, max_heating_power, 0, Ta, T_HP, compressor_speed)
                        PCM_H_Energy, left_over_thermal = charge_PCM_simple(option, PCM_H_Energy, extra_thermal_power)
                        Decision = 1
                        Scenario = 1
                    else
                        Decision = 0
                        Scenario = 0
                    end
                elseif urgency_c >= 0.7 # cooling is urgent
                    option = "C"
                    favorable = check_favorable_COP(option, input_df, PCM_C_Energy)
                    if favorable === true # charge PCM H
                        SOC, extra_thermal_power, owed_thermal_power = discharge_battery_to_HP('N', option, SOC, max_cooling_power, 0, Ta, T_HP, compressor_speed)
                        PCM_C_Energy, left_over_thermal = charge_PCM_simple(option, PCM_C_Energy, extra_thermal_power)
                        Decision = 2
                        Scenario = 5
                    else
                        Decision = 0
                        Scenario = 0
                    end
                else # urgent PCM does not exist
                    Decision = 0
                    Scenario = 0
                end
            end
        end
    end


    if PVGeneration > 0 # Day Time
        if PVGeneration < E_load
            E_load_remain = E_load - PVGeneration # need to discharge battery to meet the remaining electrical load
            SOC, loss_of_load = discharge_battery(SOC, E_load_remain) # Discharge battery to meet electrical load first
            if loss_of_load > 0 # Loss of load scenario
                loss_of_load_E = loss_of_load_E + loss_of_load
                if heating_0 == 1
                    option = "H"
                    PCM_H_Energy, owed_power = discharge_PCM(option, PCM_H_Energy, heating_power)
                    if owed_power > 0
                        loss_of_load_Q = owed_power
                        Decision = 3
                        Scenario = 4
                    else
                        Decision = 3
                        Scenario = 2
                    end
                elseif cooling_0 == 1
                    option = "C"
                    PCM_C_Energy, owed_power = discharge_PCM(option, PCM_C_Energy, cooling_power)
                    if owed_power > 0
                        loss_of_load_Q = owed_power
                        Decision = 3
                        Scenario = 8
                    else
                        Decision = 3
                        Scenario = 6
                    end
                else
                    Decision = 0
                    Scenario = 0
                end
            else # PV + battery discharge already met electrical load, PV depleted
                if heating_0 == 1 
                    option = "H"
                    SOC, extra_thermal_power, owed_thermal_power = discharge_battery_to_HP(option, option, SOC, HP_power_H, 0, Ta, T_HP, compressor_speed) # discharge battery to use HP
                    if owed_thermal_power > 0 # HP can't meet all the thermal power demand
                        PCM_H_Energy, owed_thermal_power = discharge_PCM(option, PCM_H_Energy, owed_thermal_power) # discharge PCM (controlled flow) to meet remaining thermal power
                        if owed_thermal_power > 0
                            loss_of_load_Q = owed_thermal_power
                            Decision = 3
                            Scenario = 4
                        else
                            Decision = 3
                            Scenario = 4
                        end
                    else
                        Decision = 0
                        Scenario = 9
                    end
                elseif cooling_0 == 1 # need optimization
                    option = "C"
                    SOC, extra_thermal_power, owed_thermal_power = discharge_battery_to_HP(option, option, SOC, HP_power_C, 0, Ta, T_HP, compressor_speed) # discharge battery to use HP
                    if owed_thermal_power > 0 # HP can't meet all the thermal power demand
                        PCM_C_Energy, owed_thermal_power = discharge_PCM(option, PCM_C_Energy, owed_thermal_power) # discharge PCM (controlled flow) to meet remaining thermal power
                        if owed_thermal_power > 0
                            loss_of_load_Q = owed_thermal_power
                            Decision = 3
                            Scenario = 8
                        else
                            Decision = 3
                            Scenario = 8
                        end
                    else
                        Decision = 0
                        Scenario = 10
                    end
                else # no heating or cooling Call
                    Decision = 0
                    Scenario = 0
                end
            end
        else # PV already met electrical load
            Available_PV = PVGeneration - E_load
            if heating_0 == 1 
                option = "H"
                if Available_PV >= HP_power_H # PV can fully meet HP and have left over power, guaranteed baseline heating mode
                    if Available_PV > max_heating_power # PV can support hybrid charging mode
                        PCM_H_Energy, left_over_power, left_over_thermal = hybrid_charge_PCM(option, PCM_H_Energy, Available_PV, Ta, T_HP, compressor_speed)
                        SOC, curtailed_power = charge_battery(SOC, left_over_power) # charge battery
                        # curtailment = curtailed_power
                        Decision = 1
                        Scenario = 3
                    else # PV can only support HP baseline heating mode, prioritize charging battery
                        output_thermal, left_over_power = use_HP(HP_power_H, option, Ta, T_HP, compressor_speed)
                        left_over_power = Available_PV - HP_power_H # [kW] remaining PV
                        SOC, curtailed_power = charge_battery(SOC, left_over_power) # charge battery
                        # curtailment = curtailed_power
                        Decision = 0
                        Scenario = 9
                    end
                else
                    SOC, extra_thermal_power, owed_thermal_power = discharge_battery_to_HP(option, option, SOC, HP_power_H, Available_PV, Ta, T_HP, compressor_speed) # discharge battery to use HP
                    # extra_thermal_power should be 0, or curtailed
                    if owed_thermal_power > 0 # heat pump unable to meet heating demand
                        PCM_H_Energy, owed_thermal_power = discharge_PCM(option, PCM_H_Energy, owed_thermal_power) # discharge PCM (controlled flow) to meet heating call
                        Decision = 3
                        Scenario = 4
                        if owed_thermal_power > 0
                            loss_of_load_Q = owed_thermal_power
                        else
                            loss_of_load_Q = 0
                        end
                    else # baseline heating
                        Decision = 0
                        Scenario = 9
                    end
                end
            elseif cooling_0 == 1 
                option = "C"
                if Available_PV >= HP_power_C # PV can fully meet HP and have left over power, guaranteed baseline cooling mode
                    if Available_PV > max_cooling_power # PV can support cooling hybrid charging mode
                        PCM_C_Energy, left_over_power, left_over_thermal = hybrid_charge_PCM(option, PCM_C_Energy, Available_PV, Ta, T_HP, compressor_speed)
                        SOC, curtailed_power = charge_battery(SOC, left_over_power) # charge battery
                        # curtailment = curtailed_power
                        Decision = 2
                        Scenario = 7
                    else # PV can only support HP baseline cooling mode, prioritize charging battery
                        output_thermal, left_over_power = use_HP(HP_power_C, option, Ta, T_HP, compressor_speed)
                        left_over_power = Available_PV - HP_power_C # [kW] remaining PV
                        SOC, curtailed_power = charge_battery(SOC, left_over_power) # charge battery
                        # curtailment = curtailed_power
                        Decision = 0
                        Scenario = 10
                    end
                else
                    SOC, extra_thermal_power, owed_thermal_power = discharge_battery_to_HP(option, option, SOC, HP_power_C, Available_PV, Ta, T_HP, compressor_speed) # discharge battery to use HP
                    # extra_thermal_power should be 0, or curtailed
                    if owed_thermal_power > 0 # heat pump unable to meet heating/cooling demand
                        PCM_C_Energy, owed_thermal_power = discharge_PCM(option, PCM_C_Energy, owed_thermal_power) # discharge PCM (controlled flow) to meet cooling call
                        Decision = 3
                        Scenario = 8
                        if owed_thermal_power > 0
                            loss_of_load_Q = owed_thermal_power
                        else
                            loss_of_load_Q = 0
                        end
                    else # baseline cooling
                        Decision = 0
                        Scenario = 10
                    end
                end    
            else # no heating or cooling call
                SOC, extra_power = charge_battery(SOC, Available_PV) # charge battery
                if extra_power > 0 # PV still available
                    urgency_H, urgency_C = Get_urgency_score(input_df, M_States)
                    efficiency_H, efficiency_C = Get_efficiency_score(Ta, T_HP, compressor_speed)
                    if urgency_H >= 0.7 # Heating PCM urgent
                        option = "H"
                        PCM_H_Energy, left_over_power, left_over_thermal = charge_PCM(option, PCM_H_Energy, extra_power, Ta, T_HP, compressor_speed)
                        # curtailment = left_over_power # [kW] curtailment
                        Decision = 1
                        Scenario = 1
                    elseif urgency_C >= 0.7 # Cooling PCM urgent
                        option = "C"
                        PCM_C_Energy, left_over_power, left_over_thermal = charge_PCM(option, PCM_C_Energy, extra_power, Ta, T_HP, compressor_speed)
                        # curtailment = left_over_power # [kW] curtailment
                        Decision = 2
                        Scenario = 5
                    else
                        if previous_decision == 1 # previously charging hot PCM
                            option = "H"
                            PCM_H_Energy, left_over_power, left_over_thermal = charge_PCM(option, PCM_H_Energy, extra_power, Ta, T_HP, compressor_speed)
                            # curtailment = left_over_power # [kW] curtailment
                            Decision = 1
                            Scenario = 1
                        elseif previous_decision == 2 # previously charging cold PCM
                            option = "C"
                            PCM_C_Energy, left_over_power, left_over_thermal = charge_PCM(option, PCM_C_Energy, extra_power, Ta, T_HP, compressor_speed)
                            # curtailment = left_over_power # [kW] curtailment
                            Decision = 2
                            Scenario = 5
                        else # previously not charging PCM
                            if efficiency_H >= efficiency_C # Heating PCM more efficient
                                option = "H"
                                PCM_H_Energy, left_over_power, left_over_thermal = charge_PCM(option, PCM_H_Energy, extra_power, Ta, T_HP, compressor_speed)
                                # curtailment = left_over_power # [kW] curtailment
                                Decision = 1
                                Scenario = 1
                            else # Cooling PCM more efficient
                                option = "C"
                                PCM_C_Energy, left_over_power, left_over_thermal = charge_PCM(option, PCM_C_Energy, extra_power, Ta, T_HP, compressor_speed)
                                # curtailment = left_over_power # [kW] curtailment
                                Decision = 2
                                Scenario = 5
                            end
                        end
                    end
                else # only charge battery, no extra PV available
                    Decision = 0
                    Scenario = 0
                end
            end
        end
    end

    J_States_new = DataFrame(
                "PV Generation (kW)" => [PVGeneration],
                "Lighting Plug Load (kW)" => [E_load],
                "Remain Battery Energy (kWh)" => [SOC_initial])
    
    slack = DataFrame(            
                "Curtailment (kW)" => [curtailment],
                "Loss of Load (kW)"  => [loss_of_load]
            )

    # The J_States_new consist of PV generation of the current iteration, the electrical load of the current iteration (data from current iteration of input_df)
    # the initial battery SOC at the start of the current iteration, and the curtailment/loss of load from the last iteration (solved using last iteration's PV, actual load, and initial battery SOC)

    return loss_of_load, J_States_new, slack, [Decision, Decision], Scenario
end    

function discharge_battery_to_HP(option1, option2, SOC, required_power, available_power, Ta, T_HP, compressor_speed)

    # HP only has 3 modes: hybrid charging, charging, or baseline heating/cooling
    # assume that for baseline heating/cooling, it will use HP_power_H/C (3kW)
    # assume that for hybrid charging, it will use max_heating/cooling_power (5kW)
    # assume that for charging, it will use max_heating/cooling_power (5kW)

    # required_power is the required electrical power of the heat pump
    # available_power is the current available PV power

    thermal_power = 0

    needed_discharge = required_power - available_power # [kW] power needed from the battery
    dischargable_power = min(needed_discharge/η, SOC/δt, InverterSize) # [kW]

    if option1 == "H"
        thermal_power = heating_power # [kW] amount of heat needed from the heating call
    elseif option1 == "C"
        thermal_power = cooling_power # [kW] amount of coolth needed from the cooling call
    elseif option1 == "N" # (need 2 options)
        thermal_power = 0
    end

    updated_SOC = SOC - dischargable_power * δt # [kWh] updated battery SOC
    
    discharged_power = dischargable_power * η
    
    total_power_for_HP = available_power + discharged_power 

    output_thermal, left_over_power = use_HP(total_power_for_HP, option2, Ta, T_HP, compressor_speed) # in this case, left_over_power should always be 0 since we won't discharge extra battery power

    extra_thermal_power = 0

    owed_thermal_power = 0

    if output_thermal >= thermal_power
        extra_thermal_power = output_thermal - thermal_power # [kW] this can be used to charge the PCM
    else
        owed_thermal_power = thermal_power - output_thermal # [kW] this needs to be discharged from the PCM
    end

    return updated_SOC, extra_thermal_power, owed_thermal_power
end

function charge_battery(SOC, available_power)

    chargable_power = min(available_power, (BatterySize-SOC)/δt, InverterSize) # [kW]
    updated_SOC = SOC + η * chargable_power * δt # [kWh]
    left_over_power = available_power - chargable_power # [kW]

    return updated_SOC, left_over_power
end

function discharge_battery(SOC, required_power)

    dischargable_power = min(required_power/η, SOC/δt, InverterSize) # [kW]
    
    updated_SOC = SOC - dischargable_power * δt # [kWh]
    
    loss_of_load = required_power - dischargable_power * η # [kW]

    discharged_power = dischargable_power * η
    return updated_SOC, loss_of_load
end

function charge_PCM(option, SOC, available_power, Ta, T_HP, compressor_speed)
    # use electricity to charge PCM
    # SOC: PCM state of charge
    # available_power: available electrical power for Heat Pump

    # when in Heating/Cooling Charge TES mode, max electrical power available is sent to the Heat Pump
    # the thermal output of the Heat Pump is used to charge the PCM storage.
    
    output_thermal, left_over_power = use_HP(available_power, option, Ta, T_HP, compressor_speed) # send all available PV to HP, return thermal output and left over electrical power
    if option == "H"
        chargable_power = min(output_thermal, (PCM_H_Size-SOC)/δt, PCM_H_Power_Capacity) # [kW]
        updated_SOC = SOC + chargable_power * δt # [kWh]
        left_over_thermal = output_thermal - chargable_power # [kW] the left_over_thermal technically stays in the system, but for now it is curtailed
    elseif option == "C"
        chargable_power = min(output_thermal, (PCM_C_Size-SOC)/δt, PCM_C_Power_Capacity) # [kW]
        updated_SOC = SOC + chargable_power * δt # [kWh]
        left_over_thermal = output_thermal - chargable_power # [kW] the left_over_thermal technically stays in the system, but for now it is curtailed
    else
        return 0, 0, 0
    end

    return updated_SOC, left_over_power, left_over_thermal
end

function charge_PCM_simple(option, SOC, available_thermal)
    # use heat/coolth to charge PCM
    # available_thermal: available thermal power to be used to charge PCM
    if option == "H"
        chargable_power = min(available_thermal, (PCM_H_Size-SOC)/δt, PCM_H_Power_Capacity) # [kW]
        updated_SOC = SOC + chargable_power * δt # [kWh]
        left_over_thermal = available_thermal - chargable_power # [kW] the left_over_thermal technically stays in the system, but for now it is curtailed
    elseif option == "C"
        chargable_power = min(available_thermal, (PCM_C_Size-SOC)/δt, PCM_C_Power_Capacity) # [kW]
        updated_SOC = SOC + chargable_power * δt # [kWh]
        left_over_thermal = available_thermal - chargable_power # [kW] the left_over_thermal technically stays in the system, but for now it is curtailed
    else
        return 0, 0, 0
    end

    return updated_SOC, left_over_thermal
end

function hybrid_charge_PCM(option, SOC, available_power, Ta, T_HP, compressor_speed)
    # SOC: PCM state of charge
    # available_power: available electrical power

    # when in Heating/Cooling Hybrid Charge mode, max electrical power available is sent to the Heat Pump
    # the thermal output of the Heat Pump is first used to meet the Heating/Cooling call of the ambient space
    # then, the rest of the thermal output is used to charge the PCM storage.

    # activate condition: if Available_PV > HP_heating_power, and battery is full
    output_thermal, left_over_power = use_HP(available_power, option, Ta, T_HP, compressor_speed) # send all available PV to HP, return thermal output and left over electrical power
    if option == "H"
        remaining_thermal_power = output_thermal - heating_power # [kW] heat
        chargable_power = min(remaining_thermal_power, (PCM_H_Size-SOC)/δt, PCM_H_Power_Capacity) # [kW]
        updated_SOC = SOC + chargable_power * δt # [kWh]
        left_over_thermal = remaining_thermal_power - chargable_power # [kW] the left_over_thermal technically stays in the system, but for now it is curtailed
    elseif option == "C"
        remaining_thermal_power = output_thermal - cooling_power # [kW] coolth
        chargable_power = min(remaining_thermal_power, (PCM_C_Size-SOC)/δt, PCM_C_Power_Capacity) # [kW]
        updated_SOC = SOC + chargable_power * δt # [kWh]
        left_over_thermal = remaining_thermal_power - chargable_power # [kW] the left_over_thermal technically stays in the system, but for now it is curtailed
    else
        return 0, 0, 0
    end

    return updated_SOC, left_over_power, left_over_thermal
end


function discharge_PCM(option, SOC, required_power)
    # SOC: PCM state of charge
    # required_power: required discharged thermal power
    if option == "H"
        discharged_power = min(required_power, SOC/δt, PCM_H_Power_Capacity) # [kW]
        updated_SOC = SOC - discharged_power * δt # [kWh]
        owed_power = required_power - discharged_power # [kW]
    elseif option == "C"
        discharged_power = min(required_power, SOC/δt, PCM_C_Power_Capacity) # [kW]
        updated_SOC = SOC - discharged_power * δt # [kWh]
        owed_power = required_power - discharged_power # [kW]
    else
        return 0, 0
    end

    return updated_SOC, owed_power
end


function use_HP(input_power, option, Ta, T_HP, compressor_speed)
    if option == "heating"
        COP_heating = get_COP_heating(Ta, T_HP, compressor_speed)
        actual_input_power = min(max_heating_power, input_power)
        left_over_power = input_power - actual_input_power # [kW] left over electrical power
        output_thermal = COP_heating * actual_input_power # [kW] heat gain
    elseif option == "cooling"
        COP_cooling = get_COP_cooling(Ta, T_HP, compressor_speed)
        actual_input_power = min(max_cooling_power, input_power)
        left_over_power = input_power - actual_input_power # [kW] left over electrical power
        output_thermal = - COP_cooling * actual_input_power # [kW] heat loss
    else
        return 0, 0
    end
    
    return output_thermal, left_over_power
end


function Get_urgency_score(input_df, M_States)
    
    NumTime = size(input_df)[1];
    
    Heating_Demand = 0
    Cooling_Demand = 0

    for i = 1:NumTime
        Heating_Demand = Heating_Demand + input_df.Heating_Load[i] * (1-i/NumTime)
        Cooling_Demand = Cooling_Demand + input_df.Cooling_Load[i] * (1-i/NumTime)
    end

    # Length of the M_States dataframe
    msize = size(M_States)[1];

    # Initialize the starting values with information from J_States and M_States
    PCM_C_temp = M_States[msize, :"PCM_Cold_Temp"]
    PCM_H_temp = M_States[msize, :"PCM_Hot_Temp"]

    PCM_C_Energy, PCM_H_Energy = PCM_available_energy(PCM_C_temp, PCM_H_temp)

    urgency_H = (Heating_Demand/PCM_H_Energy) / ((Heating_Demand/PCM_H_Energy) + (Cooling_Demand/PCM_C_Energy))

    urgency_C = (Cooling_Demand/PCM_C_Energy) / ((Heating_Demand/PCM_H_Energy) + (Cooling_Demand/PCM_C_Energy))

    return urgency_H, urgency_C
end    

function get_COP_heating(Ta, T_HP, compressor_speed)
    # Ta, T_HP in Celsius, compressor_speed between 0-1
    H_coeff = [10.333549711254642, 0.23633115606149574, -0.2139650595368673, -0.9691216987057145, -0.00021269966868176288, 0.0010849642966085694, 
    -0.7165709561675869, -0.0029164443083831496, -0.02624318720676654, 0.04259084351120646]
    COP_H = H_coeff[1] + H_coeff[2]*Ta + H_coeff[3]*T_HP +  H_coeff[4] * compressor_speed + H_coeff[5] * Ta^2 + H_coeff[6] * T_HP^2 + H_coeff[7] * compressor_speed^2 + 
    H_coeff[8] * (Ta * T_HP) + H_coeff[9] * (Ta * compressor_speed) + H_coeff[10] * (T_HP * compressor_speed)

    return COP_H
end

function get_COP_cooling(Ta, T_HP, compressor_speed)
    # Ta, T_HP in Celsius, compressor_speed between 0-1
    C_coeff = [0.0, -0.7246612103083289, 3.7648769496122014, -1.828441146006818, 0.008531721013710941, -0.13544346612221078, 0.6218944604023746, 
    -0.006359217720580048, 0.03675580836221801, -0.059931670614125486]
    COP_C = C_coeff[1] + C_coeff[2]*Ta + C_coeff[3]*T_HP +  C_coeff[4] * compressor_speed + C_coeff[5] * Ta^2 + C_coeff[6] * T_HP^2 + C_coeff[7] * compressor_speed^2 + 
    C_coeff[8] * (Ta * T_HP) + C_coeff[9] * (Ta * compressor_speed) + C_coeff[10] * (T_HP * compressor_speed)

    return COP_C
end

function Get_battery_SOC(PVGen, Total_Load, SOC)

    curtailment = 0
    lossofload_e = 0

    SOC = SOC* (1 - δt * BatteryLoss)

    if PVGen > Total_Load       
        SOC, curtailment = charge_battery(SOC, PVGen - Total_Load)
    else
        SOC, lossofload_e = discharge_battery(SOC, Total_Load - PVGen)
    end

    return SOC, curtailment, lossofload_e
end


function PCM_available_energy(PCM_C_temp, PCM_H_temp)

    # 10 degrees each way outside of the latent range
    
    PCM_C_Energy = 0
    PCM_H_Energy = 0

    if PCM_C_temp > 12
        PCM_C_Energy = (22 - PCM_C_temp) * 2050 * 271.14 / (3.6*10^6)
    end

    if PCM_C_temp >= 10 && PCM_C_temp <=12
        PCM_C_Energy = ((12 - PCM_C_temp) * (126000/2) + 10 * 2050) * 271.14 / (3.6*10^6)
    end

    if PCM_C_temp < 10
        PCM_C_Energy = ((10 - max(0, PCM_C_temp)) * 2050 + 126000 + 10 * 2050) * 271.14 / (3.6*10^6)
    end

    if PCM_H_temp < 47
        PCM_H_Energy = (PCM_H_temp - 37) * 3150 * 179.2 / (3.6*10^6)
    end

    if PCM_H_temp >= 47 && PCM_H_temp <=49
        PCM_H_Energy = ((PCM_H_temp - 47) * (180800/2) + 10 * 3150) * 179.2 / (3.6*10^6)
    end

    if PCM_H_temp > 49
        PCM_H_Energy = ((min(59, PCM_H_temp) - 49) * 3150 + 180800 + 10 * 3150) * 179.2 / (3.6*10^6)
    end

    PCM_H_Energy_corrected = max(0, PCM_H_Energy) # [kWh]
    PCM_C_Energy_corrected = max(0, PCM_C_Energy) # [kWh]

    return PCM_C_Energy_corrected, PCM_H_Energy_corrected
end

function Get_efficiency_score(Ta, T_HP, compressor_speed)
    # The efficiency score should only be the efficiency of the charging mode, not the hybrid charging mode.
    # Since hybrid charging mode is always decided by the heating/cooling call itself.

    # COP of HP not considered, but should be considered later.
    HP_SP_C_charge = 8 
    HP_SP_H_charge = 55 
    
    efficiency_score_H = 1/(abs(T_HP - HP_SP_H_charge)+0.1)
    efficiency_score_C = 1/(abs(T_HP - HP_SP_C_charge)+0.1)

    return efficiency_score_H, efficiency_score_C
end 

function check_favorable_COP(option, input_df, PCM_SOC)
    # first find the horizon until the solar is available again.
    # Find the first index where PV > 0
    first_positive_index = findfirst(x -> x > 0, input_df.PV)

    # Split the DataFrame
    if first_positive_index !== nothing
        df_night = input_df[1:first_positive_index-1, :]  # Before first positive value
        df_after = input_df[first_positive_index:end, :]   # From first positive value onwards
    else
        df_night = input_df  # If no positive values, all goes in df_night
        df_after = DataFrame()  # Empty DataFrame
    end

    n = nrow(df_night)
    if n == 1
        return false
    end

    current_Ta = df_night[1, :Ta] 

    # how much should the battery be discharging to use the better COP? what about forecast uncertainty? should we throw everything we have assuming that the 
    # future heating demand is 100% there and the COP 100% gets worse?
    if option == "H"

        # check Ta now, and all the Ta until the sun comes back up.
        # if there this is the highest Ta, and there is heating load before the sun comes up, return true
        # if there this is the highest Ta until another higher Ta before the sun comes up, and there is heating load in between, return true
        # else return false

        # sum up all heating demand of the night
        # check current SOC of PCM H
        # if PCM H SOC <= 1.2 * all heating demand
        # if the current Ta is absolutely the highest,
        # charge, and then update PCM H SOC

        # if the current Ta is the highest until a later timestep,
        # cut the horizon and consider the new horizon
        # if PCM H SOC <= 1.2 * heating demand in the new horizon
        # charge, and then update PCM H SOC

        index = findfirst(>(current_Ta), df_night.Ta)
        
        if index === nothing
            index = n
        end
        
        if index == 2
            return false
        end

        next_horizon = df_night[1:index-1, :]

        future_heating_load = sum(df_night[1:index-1, :Heating_Load])*δt
        if PCM_SOC <= future_heating_load * 1.2
            return true
        else
            return false
        end
    elseif option == "C"
        index = findfirst(<(current_Ta), df_night.Ta)
        
        if index === nothing
            index = n
        end
        
        if index == 2
            return false
        end

        next_horizon = df_night[1:index-1, :]

        future_cooling_load = sum(df_night[1:index-1, :Cooling_Load])*δt
        if PCM_SOC <= future_cooling_load * 1.2
            return true
        else
            return false
        end
    end
end


#=
function optimize_heat_pump_schedule_H(df::DataFrame, SOC_0, PCM_H_Energy, T_HP, compressor_speed)
    
    # the assumption is that we can control exactly how much we could discharge the battery to use the HP at a given point.
    # therefore, H2HP or H2C can be optimized.

    n = nrow(df)
    model = Model(Gurobi.Optimizer)  # Use GLPK solver (or replace with Gurobi, CPLEX, etc.)

    # Decision variables
    @variable(model, H2HP[1:n] >= 0)  # [kW] Electric power used by heat pump for heating
    @variable(model, PCM_H_SOC[1:n] >= 0) # [kWh] Heat stored in thermal storage at each time step
    @variable(model, PCM_H2H[1:n] >= 0) # [kW] Heat discharged from storage at each time step
    @variable(model, SOC[1:n] >= 0) # [kW] Battery SOC

    # Objective: Minimize total electricity usage
    @objective(model, Min, sum(H2HP[t] for t in 1:n))

    @expression(model, COP_H[t], get_COP_heating(df[t, :Ta], T_HP, compressor_speed))
    # Constraints
    for t in 1:n
        # Heat balance: Heating load must be met by HP or storage discharge
        @constraint(model, H2HP[t] * COP_H[t] + PCM_H2H[t] >= df[t, :Heating_Load])

        # Thermal storage balance
        if t == 1
            # Initial storage is zero (or you can set it to a value)
            @constraint(model, PCM_H_SOC[t] == PCM_H_Energy)
        else
            # Storage at time t = storage at time t-1 + heat generated - heat discharged
            @constraint(model, PCM_H_SOC[t] == PCM_H_SOC[t-1] + H2HP[t] * COP_H[t] - PCM_H2H[t])
        end

        # Storage capacity constraint
        @constraint(model, PCM_H_SOC[t] <= PCM_H_Size)

        # Heat pump capacity constraint
        @constraint(model, H2HP[t] <= max_heating_power)

        # battery power capacity constraint (should automatically satisfy)
        # @constraint(model, H2HP[t] <= InverterSize)
        # Thermal storage balance
        if t == 1
            # Initial storage is zero (or you can set it to a value)
            @constraint(model, SOC[t] == SOC_0)
        else
            # Storage at time t = storage at time t-1 + heat generated - heat discharged
            @constraint(model, SOC[t] == SOC[t-1] - (H2HP[t]/η) * δt)
        end
    end

    # Solve the model
    optimize!(model)

    # Extract results
    HP_Action = value.(H2HP[1])
    Updated_PCM_H_SOC = value.(PCM_H_SOC[1])
    PCM_Discharge = value.(PCM_H2H[1])

    return HP_Action, PCM_Discharge, Updated_PCM_H_SOC
end

function optimize_heat_pump_schedule_C(df::DataFrame, PCM_C_Energy, T_HP, compressor_speed)
    
    n = nrow(df)
    model = Model(Gurobi.Optimizer)  # Use GLPK solver (or replace with Gurobi, CPLEX, etc.)

    # Decision variables
    @variable(model, H2C[1:n] >= 0)  # [kW] Electric power used by heat pump for heating
    @variable(model, PCM_C_SOC[1:n] >= 0) # [kWh] Heat stored in thermal storage at each time step
    @variable(model, PCM_C2H[1:n] >= 0) # [kW] Heat discharged from storage at each time step

    # Objective: Minimize total electricity usage
    @objective(model, Min, sum(H2C[t] for t in 1:n))

    @expression(model, COP_C[t], get_COP_cooling(df[t, :Ta], T_HP, compressor_speed))
    # Constraints
    for t in 1:n
        # Heat balance: Heating load must be met by HP or storage discharge
        @constraint(model, H2C[t] * COP_C[t] + PCM_C2H[t] >= df[t, :Cooling_Load])

        # Thermal storage balance
        if t == 1
            # Initial storage is zero (or you can set it to a value)
            @constraint(model, PCM_C_SOC[t] == PCM_C_Energy)
        else
            # Storage at time t = storage at time t-1 + heat generated - heat discharged
            @constraint(model, PCM_C_SOC[t] == PCM_C_SOC[t-1] + H2C[t] * COP_C[t] - PCM_C2H[t])
        end

        # Storage capacity constraint
        @constraint(model, PCM_C_SOC[t] <= PCM_C_Size)

        # Heat pump capacity constraint
        @constraint(model, H2C[t] <= max_cooling_power)
    end

    # Solve the model
    optimize!(model)

    # Extract results
    HP_Action = value.(H2C[1])
    Updated_PCM_C_SOC = value.(PCM_C_SOC[1])
    PCM_Discharge = value.(PCM_C2H[1])

    return HP_Action, PCM_Discharge, Updated_PCM_C_SOC
end

function optimize_pcm_h_charging(night_df::DataFrame, PCM_H_SOC::Float64)
    # 1. Define the night horizon (assume night_df is already the PV-off to PV-on period)
    # If not pre-defined, you would filter the full dataset to where PV generation is zero.
    # For example (if PV data available): night_df = df[df.PV .== 0, :] 
    # Here, we assume night_df is the relevant slice.
    local_horizon = night_df  # working copy of the horizon DataFrame
    
    # 2. Compute total heating demand over the current horizon
    total_demand = sum(local_horizon.Heating_Load)
    
    # Check if storage is more than sufficient for the horizon
    if PCM_H_SOC > 1.2 * total_demand
        # No need to charge at all during this horizon
        return PCM_H_SOC, false, local_horizon  # storage unchanged, no charge, horizon unchanged
    end
    
    # 3 & 4. Decide on charging with potential horizon adjustments
    charge_now = false
    while true
        # Recalculate demand for current horizon segment (in case horizon was adjusted)
        total_demand = sum(local_horizon.Heating_Load)
        
        # If even for the adjusted horizon the storage is sufficient, no charge needed for now
        if PCM_H_SOC > 1.2 * total_demand
            # We can cover this horizon with existing storage; no immediate charge
            break  # exit loop with charge_now = false
        end
        
        # Determine ambient temperature conditions
        current_Ta = local_horizon.Ta[1]                   # ambient temp at start of horizon
        max_Ta_horizon = maximum(local_horizon.Ta)         # peak ambient temp in this horizon
        
        if current_Ta >= max_Ta_horizon
            # 4a. Current ambient is the highest in the horizon -> charge now
            charge_now = true
            # 5. Compute how much to charge (to reach 1.2 × horizon demand)
            local_needed = 1.2 * total_demand - PCM_H_SOC
            if local_needed < 0 
                local_needed = 0 
            end
            PCM_H_SOC += local_needed
            # (Optionally cap PCM_H_SOC to its max capacity if known, e.g., PCM_H_SOC = min(PCM_H_SOC, PCM_capacity))
            break  # charging done, exit loop
        else
            # 4b. A higher ambient temperature occurs later -> adjust horizon to wait until then
            # Find the next time index where ambient temperature exceeds the current_Ta
            next_higher_idx = findfirst(x -> x > current_Ta, local_horizon.Ta)
            # findfirst will return the first index (including the first row) that meets the condition.
            # We need the first index AFTER the start, so ensure it's not index 1:
            if next_higher_idx == 1
                next_higher_idx = findnext(x -> x > current_Ta, local_horizon.Ta, 2)
            end
            
            if next_higher_idx === nothing
                # No higher temperature found (edge case) – treat current as highest (should not happen here)
                charge_now = true
                local_needed = 1.2 * total_demand - PCM_H_SOC
                if local_needed < 0 
                    local_needed = 0 
                end
                PCM_H_SOC += local_needed
                break
            end
            
            # Redefine the horizon to end just before the next higher-Ta time
            local_horizon = local_horizon[1:next_higher_idx-1, :]
            # Loop will repeat: recalc demand for this shorter horizon and check again
            continue
        end
    end
    
    # 6. Return the updated SOC, the charging decision, and the horizon used for decision
    return PCM_H_SOC, charge_now, local_horizon
end

function Get_battery_SOC_old(PVGen_0, Total_Load_0, B_SOC_0)

    # change rulebased

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
    @constraint(m_0, InStorageBattery_past[2] == InStorageBattery_past[1] * (1 - δt * BatteryLoss) + (PV2B_past * η - B2H_past)*δt); # [kWh]
    @constraint(m_0, Total_Load_0 == PV2H_past * η_PVIV + B2H_past * η) # [kW]

    # Battery discharging constraint, node at battery 
    @constraint(m_0, InStorageBattery_past[1] == B_SOC_0); # [kWh]

    # Battery discharging constraint, node at battery 
    @constraint(m_0, B2H_past * δt <= InStorageBattery_past[1]); # [kWh]
    
    # Battery power inverter constraint, node at battery (inverter power constraint) 
    @constraint(m_0, (B2H_past + PV2B_past) <= InverterSize); # [kW]
    
    # Battery storage size constraint, node at battery
    @constraint(m_0, [t=1:2], InStorageBattery_past[t] <= BatterySize); # [kWh]
    
    # Battery storage max discharge constraint, node at battery, always at least 20% full 
    @constraint(m_0, [t=1:2], InStorageBattery_past[t] >= BatterySize * (1-MaxDischarge)); # [kWh]

    optimize!(m_0);
    
    # Actual Battery SOC to start with
    InStorageBattery_1 = value.(InStorageBattery_past[2]) # [kWh]
    return InStorageBattery_1
end

=#