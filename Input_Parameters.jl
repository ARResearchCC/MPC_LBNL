# Input Parameters 
urgency_threshold = 0.4 # 0.3 and 0.7

Initial_Ti = 25.5 # [°C]
Initial_PV_Gen = 0 # [kW]
Initial_P_0 = 0 # [kW]
Initial_B_SOC = 1 # full battery to start with
Initial_Curtailment = 0 # [kW]
Initial_Loss_of_Load = 0 # [kW]

Inital_PCM_H_Temp = 48 # [°C]
Inital_PCM_C_Temp = 11 # [°C]

Intial_PCM_C_SOC = 0.5
Intial_PCM_H_SOC = 0.5

# Interior Standard
begin
    # Temperature standard
    # Temperature setpoints control the direct heating or cooling of the space through heat pump or PCM discharge.
    zone_temp_cooling_setpoint = 28 # [°C]
    zone_temp_heating_setpoint = 23 # [°C]
    zone_temp_setpoint_delta = 2 # [°C]

    # Ventilation
    Ventilation = 15 # [CFM/PPL] Building Standard
    Ra = 0.06 # [CFM/ft^2] Area Outdoor Air Rate; Source:CEE226E Week 5 - Energy Modeling Questions - Slide 9
    Rp = 5 # [CFM/PPL] People Outdoor Air Rate; Source:CEE226E Week 5 - Energy Modeling Questions - Slide 9   
end

# Lighting, Plugs, and Occupancy
begin
    PeakLighting = 0.1 # [KW] Dan's suggestion
    PeakPlugLoad = 0.1 # [KW] Computer = 40W, Phone = 10W, 2*(40 + 10) = 100 [W] = 0.1 [kW]
    MaxOccupancy = 4 # [PPL] 4-MAN Office Room Plan
    PersonLatentHeat = 200; # [BTU/hr/PPL]  CEE226E Slide
    PersonSensibleHeat = 300; # [BTU/hr/PPL]  CEE226E Slide
    TotalPersonHeat = PersonSensibleHeat + PersonLatentHeat # [BTU/hr/PPL]
end

# Device Capacity Parameters
PVSize = 27; # [kW] PV DC Power Capacity
BatterySize = 15; # [kWh] Battery Energy Capacity
InverterSize = 15 # [kW] Max Continuous AC Output Power

# Standard Operating Power of Heat Pump
HP_power_H = 3 # [kW] default constant electrical power consumption for heat pump (heating)
HP_power_C = 3 # [kW] default constant electrical power consumption for heat pump (cooling)

heating_power = 3 # [kW] default heating call thermal power requirement
cooling_power = 3 # [kW] default cooling call thermal power requirement

max_heating_power = 5 # [kW] max heat pump heating capacity
max_cooling_power = 5 # [kW] max heat pump cooling capacity

# Full Capacity
PCM_H_Size_Full = 12.135822222222222 # [kWh]
PCM_C_Size_Full = 12.577883333333334 # [kWh]

# Nominal Capacity
PCM_H_Size = 9 # [kWh]
PCM_C_Size = 9 # [kWh]

PCM_H_Power_Capacity = 5 # [kW]
PCM_C_Power_Capacity = 5 # [kW]

# Battery parameters
BatteryLoss = 0.00001 # [/hr] Battery Leakage Rate ??
MaxDischarge = 0.8 # [1] Max Depth of Discharge in Battery (80%)
η = 0.98 # [1] Battery Inverter Efficiency

# Standard Operating Power of Control System, Air Handling Unit (AHU), and Pump
P_Controls = 0.02 # [kW]
P_AHU = 0 # [kW]
P_Pumps = 0.04 # [kW]

# Solar PV Parameters
noct_installed = 45 # [°C] The “installed” nominal operating cell temperature. PVWatts assumes this value to be 45 C for rack-mounted arrays and 49 C for roof mount systems with restricted air flow around the module.
module_height = 5 # [m] The height above ground of the center of the module. The PVWatts default is 5.0.
module_width = 0.31579 # [m] Module width. The default value of 0.31579 meters in combination with the default module_length gives a hydraulic diameter of 0.5.
module_length = 1.2 # [m] Module length. The default value of 1.2 meters in combination with the default module_width gives a hydraulic diameter of 0.5.
module_emissivity = 0.84 # [1] The effectiveness of the module at radiating thermal energy.
module_absorption = 0.83 # [1] The fraction of incident irradiance that is converted to thermal energy in the module. 
module_surface_tilt = 27 # [DEGREES] Module tilt from horizontal. If not provided, the default value of 30 degrees is used.
Γ_t = -0.47/100 # [1/°C] Temperature coefficient for standard module
T_ref = 25 # [°C] Reference cell temperature
η_PV = 0.86 # [1] PV DC efficiency after system loss
P_dc0 = 1 # [kW/kW] Rated capacity at standard conditions
η_PVIV = 0.94 # [1] PV(DC) to Home(AC) inverter efficiency

# Equipments Parameters
begin       
    # Heat Pump heating and cooling coefficients of performance
    
    # HP_a = 6.08
    # HP_b = -0.0941
    # HP_c = 0.000464
    
    COP_H = 3.7 # COP of heating (currently constant)
    COP_C = 3.7 # COP of cooling (currently constant)

    # C_HP_OP = 0.02 * Cap_HP_H # [$/(kW*YR)] Operational Cost of Heat Pump

    # η_PCM_H = 1/24 # [/hr] PCM Heating Storage Leakage Rate
    # η_PCM_C = 1/24 # [/hr] PCM Cooling Storage Leakage Rate

    C_PCM_H_OP = 0.02 * PCM_H_Size # [$/(kWh*YR)] Operational Cost of PCM Heating Storage
    C_PCM_C_OP = 0.02 * PCM_C_Size # [$/(kWh*YR)] Operational Cost of PCM Cooling Storage

    PCM_H_discharge_rate = 2 # [kW] default constant heat discharging rate of PCM Heating Storage
    PCM_C_discharge_rate = 2 # [kW] default constant heat discharging rate of PCM Cooling Storage
end

############ Declare Parameters ############
begin
    # Berg Envelope Parameters
    begin
        # unit modification has been made
        L_wall = 6.1 # [m] Length of the Wall
        H_wall = 2.4 # [m] Height of the Wall
        R_wall = 2.6 # [m^2·K/W] Thermal Resistance of the Wall  It's based on time units of seconds
        R_floor = 3.0 # [m^2·K/W] Thermal Resistance of the Floor
        H_Ceiling = 2.0 # [m] Interior Height of the Ceiling (78.5 inch) 
        Area = L_wall * L_wall # [m^2] Floor/Ceiling Area 
        UA = L_wall*H_wall*4*(1/R_wall) + Area*2*(1/R_floor) # [W/K]
        Volume = L_wall * L_wall * H_Ceiling # [m^3]
        e_Berg = 0.75 # fibreglass https://www.thermoworks.com/emissivity-table/
        Berg_tilt = 15
    end
    # Geographical Parameters
    begin
        standard_meridian_longitude = -120; # [DEGREES] Longitude of standard meridian
        longitude = -122.4286 # [DEGREES] Longitude of site
        latitude = 37.4636 # [DEGREES] Latitude of site
    end
    # Infiltration Parameters
    begin
        # Stack coefficient for building(house height = 1 story)
        Cs = 0.000145 # [(L/s)^2/(cm^4*°K)]

        # Wind coefficient for building(house height =1 story, shelter class=1(no obstructions or local shielding))
        Cw = 0.000319 # [(L/s)^2/(cm^4*(m/s)^2)] (between 0.000319 and 0.000246 (house height =1 story, shelter class=2(typical shelter for an insolated rural house)))

        # Effective leakage area measured during Blower Door Test 
        # ELA =  47.1 # [in^2] (between 38.3(Dan's interpolation) and 47.1(Jessie's interpolation))
        Al =  303.9 # [cm^2]
        T_d = 273.15 - 1.67 # [°K] 99% Design Temperature for Half Moon Bay 
    end
    # Daylight Saving
    begin
        
        dls_start = DateTime(2024, 3, 10, 2)
        dls_end = DateTime(2024, 11, 3, 0)
        timezone = "America/Los_Angeles"  # Specify local timezone
        
    end
end
# Environmental Parameters
begin
    wind_height = 9.144 # [m] The height above ground at which wind speed is measured. The PVWatts default is 9.144 m.
    albedo = 0.18 # [1] albedo of grass
    n_air = 1 # [1] Index of reflaction of air
    n_AR = 1.3 # [1] Index of reflaction of AR coating
    n_glass = 1.526 # [1] Index of reflaction of glass as in PVWatt model for standard modules
end
# Infiltration Parameters
begin
    # Stack coefficient for building(house height = 1 story)
    Cs = 0.000145 # [(L/s)^2/(cm^4*K)]
    # Wind coefficient for building(house height =1 story, shelter class=1(no obstructions or local shielding))
    Cw = 0.000319 # [(L/s)^2/(cm^4*(m/s)^2)] (between 0.000319 and 0.000246 (house height =1 story, shelter class=2(typical shelter for an insolated rural house)))
    
    # Effective leakage area measured during Blower Door Test 
    ELA =  38.3 # [in^2] (between 38.3(Dan's interpolation) and 47.1(Jessie's interpolation))
    Al = ELA * 6.4516 # [cm^2]

    T_d = 273.15 - 1.67 # [°K] 99% Design Temperature for Half Moon Bay
    T_indoor_constant = 22 # [°C] Constant Indoor Temperature (for the simplicity of a linear model)
end