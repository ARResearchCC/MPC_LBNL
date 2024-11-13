# Generate Schedules Fred Fan 10/28/2024
############ Initialize Tools ############
begin
    import Pkg;
    # Initialize JuMP to allow mathematical programming models
    # Add Packages if you are running this for the first time
    
    # Add required packages
    # Pkg.add(["JuMP", "CSV", "DataFrames", "Dates", "XLSX", "StatsBase", "Statistics"])

    # Import necessary libraries
    using JuMP, CSV, DataFrames, Dates, XLSX, StatsBase, Statistics
    include("Input_Parameters.jl")
end

############ Generate Schedules ############

# Instructions
begin
    # The are two schedules, format: Schedules[Lighting, Plugs, Occupancy]
    # All variables are unitless(0.5 = 50%). They are then used to multiply peak lighting/occupancy/plugs heat gain in the code.

    # SimpleSchedule:
    # About lighting schedule: full lighting intensity from 6 am to 10 pm daily
    # About plugs schedule: full plugs intensity from 6 am to 10 pm daily
    # About occupancy schedule: full occupancy from 6 am to 10 pm daily

    # ComplexSchedule:
    # About lighting schedule: full lighting intensity from 6 am to 10 pm daily
    # About plugs schedule: equipments (PC and phone) are charged from 7 pm to 12 am daily (5 hours)
    # About occupancy schedule: 50% occupany from 8 am to 8 pm, 100% occupancy during night time.
end

function generate_schedules(f_run::Int, total_intervals, schedule)
    if schedule == "simple"
        return simpleschedule(f_run, total_intervals)
    elseif schedule == "complex"
        return complexschedule(f_run, total_intervals)
    end
end

function simpleschedule(f_run::Int, total_intervals)
    # Initialize schedules
    SimpleSchedule = zeros(total_intervals, 3)
    # Fill SimpleSchedule
    for i in 0:364
        for j in 10*f_run:(17*f_run-1)
            SimpleSchedule[i*24*f_run + j + 1, :] .= 1
        end
    end
    schedule = DataFrame()
    schedule.Lighting = SimpleSchedule[:, 1] * PeakLighting # [kW]
    schedule.Plugs = SimpleSchedule[:, 2] * PeakPlugLoad # [kW]
    schedule.Occupancy = SimpleSchedule[:, 3] * (MaxOccupancy * TotalPersonHeat/3412.14) # [kW]
    return schedule
end


function complexschedule(f_run::Int, total_intervals)
    # Initialize schedules
    ComplexSchedule = zeros(total_intervals, 3)
    
    # Fill ComplexSchedule
    for i in 0:364
        for j in 1*f_run:(6*f_run-1)
            ComplexSchedule[i*24*f_run + j + 1, 1:2] .= 0
            ComplexSchedule[i*24*f_run + j + 1, 3] = 1
        end
        for j in 7*f_run:(8*f_run-1)
            ComplexSchedule[i*24*f_run + j + 1, 1] = 1
            ComplexSchedule[i*24*f_run + j + 1, 2] = 0
            ComplexSchedule[i*24*f_run + j + 1, 3] = 1
        end
        for j in 9*f_run:(19*f_run-1)
            ComplexSchedule[i*24*f_run + j + 1, 1] = 1
            ComplexSchedule[i*24*f_run + j + 1, 2] = 0
            ComplexSchedule[i*24*f_run + j + 1, 3] = 0.5
        end
        for j in 20*f_run:(20*f_run+f_run-1)
            ComplexSchedule[i*24*f_run + j + 1, 1:2] .= 1
            ComplexSchedule[i*24*f_run + j + 1, 3] = 0.5
        end
        for j in 21*f_run:(22*f_run-1)
            ComplexSchedule[i*24*f_run + j + 1, :] .= 1
        end
        for j in 23*f_run:(24*f_run-1)
            ComplexSchedule[i*24*f_run + j + 1, 1] = 0
            ComplexSchedule[i*24*f_run + j + 1, 2] = 1
            ComplexSchedule[i*24*f_run + j + 1, 3] = 1
        end
    end
    schedule = DataFrame()
    schedule.Lighting = ComplexSchedule[:, 1] * PeakLighting # [kW]
    schedule.Plugs = ComplexSchedule[:, 2] * PeakPlugLoad # [kW]
    schedule.Occupancy = ComplexSchedule[:, 3] * (MaxOccupancy * TotalPersonHeat/3412.14) # [kW]
    return schedule
end
