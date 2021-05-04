"""
    ModelParamSet(velocity_sampling, farm_model_with_ti, farm_model_no_ti, turbine_design, turbine_op, turbine_power_models, turbine_ct_models, wind_resource_model)

Container holding all parameters for the wind farm analysis model.

# Arguments
- `velocity_sampling::RotorSampling`: Container holding rotor sampling points
- `farm_model_with_ti::WindFarmModelSet`: FlowFarm container for objects defining models to use in wind farm calculations (with turbulence intensity)
- `farm_model_no_ti::WindFarmModelSet`: FlowFarm container for objects defining models to use in wind farm calculations (without turbulence intensity)
- `turbine_design::TurbineType`: Container holding design parameters for the turbines
- `turbine_op`: Container holding parameters for the operating state of the turbines
- `turbine_power_models::AbstractPowerModel`: FlowFarm container holding the power models for each turbine
- `turbine_ct_models::AbstractThrustCoefficientModel`: Flowfarm container holding the thrust coefficient models for each turbine
- `wind_resource_model::AbstractWindResourceModel`: Flowfarm container holding the wind resource model
"""
struct ModelParamSet{VS, FMWTI, FMNTI, TD, TO, TPM, TCM, WRM}
    # farm parameters and models
    velocity_sampling::VS
    farm_model_with_ti::FMWTI
    farm_model_no_ti::FMNTI

    # turbine parameters and models
    turbine_design::TD
    turbine_op::TO
    turbine_power_models::TPM
    turbine_ct_models::TCM
    
    # wind resource model
    wind_resource_model::WRM
end

"""
    model_set(_wake_model, _turbine_type, nturbines, _windrose, _ndirs, _nspeeds)

Creates a ModelParamSet container that holds all parameters for the wind farm analysis model.

# Arguments
- `_wake_model::String`: name of the wake model
- `_turbine_type::String`: name of turbine type
- `nturbines::Int64`: number of turbines in the farm
- `_windrose::String`: name of wind rose
- `_ndirs::Int64`: number of directions in the wind rose
- `_nspeeds::Int64`: number of speeds in the wind rose
"""
function model_set(_wake_model, _turbine_type, nturbines, _windrose, _ndirs, _nspeeds)
        
    #################################################################################
    # SET FARM MODELS AND PARAMETERS
    #################################################################################

    # set rotor sample points
    velocity_sampling = RotorSampling()

    # wake deficit model
    eval(Meta.parse("wake_deficit_model = ff.$(_wake_model)()"))
    # wake_deficit_model = ff.GaussYawVariableSpread()

    # wake deflection model
    wake_deflection_model = ff.GaussYawDeflection()

    # wake combination model
    wake_combination_model = ff.LinearLocalVelocitySuperposition()

    # initialize model set
    farm_model_no_ti = ff.WindFarmModelSet(wake_deficit_model, wake_deflection_model, wake_combination_model, ff.LocalTIModelNoLocalTI())
    farm_model_with_ti = ff.WindFarmModelSet(wake_deficit_model, wake_deflection_model, wake_combination_model, ff.LocalTIModelMaxTI())


    #################################################################################
    # SET TURBINE MODELS AND PARAMETERS
    #################################################################################

    # turbine type
    eval(Meta.parse("turbine_design = $(_turbine_type)($(nturbines))"))
    # turbine_design = VestasV80_2MW(layout.nturbines)

    # operating condition
    turbine_op = NoYaw(nturbines)

    # get turbine power and thrust models
    turbine_power_model, turbine_ct_model = get_turbine_power_thrust_models(turbine_design)


    #################################################################################
    # SET WIND RESOURCE MODEL
    #################################################################################
    eval(Meta.parse("windrose = $(_windrose)($(_ndirs), $(_nspeeds))"))
    # windrose = NantucketWindRose()
    wind_resource_model = get_wind_resource_model(windrose)


    #################################################################################
    # SET PARAMETERS
    #################################################################################

    model_param = ModelParamSet(velocity_sampling, farm_model_with_ti, farm_model_no_ti, turbine_design, turbine_op,
                    turbine_power_model, turbine_ct_model, wind_resource_model)

    return model_param
end