abstract type TurbineType end

"""
    VestasV80_2MW(hub_height, rotor_diameter, cut_in_speed, cut_out_speed, rated_speed, rated_power, generator_efficiency)

Container holding design parameters for the Vestas V80 2MW turbine

# Arguments
- `hub_height::Float64`
- `rotor_diameter::Float64`
- `cut_in_speed::Float64`
- `cut_in_speed::Float64`
- `rated_speed::Float64`
- `rated_power::Float64`
- `generator_efficiency::Float64`
"""
struct VestasV80_2MW{F} <: TurbineType
    hub_height::F
    rotor_diameter::F
    cut_in_speed::F
    cut_out_speed::F
    rated_speed::F
    rated_power::F
    generator_efficiency::F
end
VestasV80_2MW(n::Int64) = VestasV80_2MW(70.0 .+ zeros(n), 80.0 .+ zeros(n), 4.0 .+ zeros(n), 25.0 .+ zeros(n), 16.0 .+ zeros(n), 2.0e6 .+ zeros(n), 0.944 .+ zeros(n))  # default: Vestas V80 2MW turbine


"""
    NREL_5MW(hub_height, rotor_diameter, cut_in_speed, cut_out_speed, rated_speed, rated_power, generator_efficiency)

Container holding design parameters for the Vestas V80 2MW turbine

# Arguments
- `hub_height::Float64`
- `rotor_diameter::Float64`
- `cut_in_speed::Float64`
- `cut_in_speed::Float64`
- `rated_speed::Float64`
- `rated_power::Float64`
- `generator_efficiency::Float64`
"""
struct NREL_5MW{F} <: TurbineType
    hub_height::F
    rotor_diameter::F
    cut_in_speed::F
    cut_out_speed::F
    rated_speed::F
    rated_power::F
    generator_efficiency::F
end
NREL_5MW(n::Int64) = NREL_5MW(90.0 .+ zeros(n), 126.4 .+ zeros(n), 3.0 .+ zeros(n), 25.0 .+ zeros(n), 11.4 .+ zeros(n), 5.0e6 .+ zeros(n), 0.944 .+ zeros(n))


"""
    TurbineOperatingParams(yaw)

Container holding parameters for the turbine operating condition

# Arguments
- `yaw::Float64`: yaw angle of the turbine w.r.t. free-stream wind direction
"""
struct NoYaw{AF}
    yaw::AF
end
NoYaw(n::Int64) = NoYaw(0.0 .+ zeros(n))


"""
    get_turbine_power_thrust_models(turbine_design::VestasV80_2MW)

Creates power and thrust models (used by FlowFarm) for the Vestas V80 2MW turbine.

# Arguments
- `turbine_design::VestasV80_2MW`: Container holding design parameters for the Vestas V80 2MW turbine
"""
function get_turbine_power_thrust_models(turbine_design::VestasV80_2MW)

    # get number of turbines
    nturbines = length(turbine_design.hub_height)

    # get turbine data
    turbinepowerdata = readdlm("../data/input-files/vestas_v80_2MW_power_curve.txt", skipstart=1)
    velpoints_power = turbinepowerdata[:,1]
    powerpoints = turbinepowerdata[:,2] * 1e6
    turbinethrustdata = readdlm("../data/input-files/vestas_v80_2MW_thrust_coefficient_curve.txt", skipstart=1)
    velpoints_thrust = turbinethrustdata[:,1]
    ctpoints = turbinethrustdata[:,2]

    # turbine power model
    turbine_power_model_single = ff.PowerModelPowerPoints(velpoints_power, powerpoints) 
    turbine_power_model = Vector{typeof(turbine_power_model_single)}(undef, nturbines)
    for i = 1:nturbines
        turbine_power_model[i] = turbine_power_model_single
    end

    # turbine thrust model
    turbine_ct_model_single = ff.ThrustModelCtPoints(velpoints_thrust, ctpoints)   
    turbine_ct_model = Vector{typeof(turbine_ct_model_single)}(undef, nturbines)
    for i = 1:nturbines
        turbine_ct_model[i] = turbine_ct_model_single
    end

    return turbine_power_model, turbine_ct_model
end

"""
    get_turbine_power_thrust_models(turbine_design::NREL_5MW)

Creates power and thrust models (used by FlowFarm) for the NREL 5MW reference turbine.

# Arguments
- `turbine_design::NREL_5MW`: Container holding design parameters for the NREL 5MW reference turbine
"""
function get_turbine_power_thrust_models(turbine_design::NREL_5MW)

    # get number of turbines
    nturbines = length(turbine_design.hub_height)

    # get turbine data
    turbinedata = readdlm("../data/input-files/NREL5MWCPCT.txt", skipstart=2)
    velpoints = turbinedata[:,1]
    cppoints = turbinedata[:,2]
    ctpoints = turbinedata[:,3]

    # turbine power model
    turbine_power_model_single = ff.PowerModelCpPoints(velpoints, cppoints) 
    turbine_power_model = Vector{typeof(turbine_power_model_single)}(undef, nturbines)
    for i = 1:nturbines
        turbine_power_model[i] = turbine_power_model_single
    end

    # turbine thrust model
    turbine_ct_model_single = ff.ThrustModelCtPoints(velpoints, ctpoints)   
    turbine_ct_model = Vector{typeof(turbine_ct_model_single)}(undef, nturbines)
    for i = 1:nturbines
        turbine_ct_model[i] = turbine_ct_model_single
    end

    return turbine_power_model, turbine_ct_model
end
