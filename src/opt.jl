"""
    optimize_farm_layout(final_layout_path, opt_info_directory, layout_param, model_param, opt_algorithm::SnoptWECAlgorithm, boundary::CircleBoundary)

Runs the turbine layout optimization for a circular boundary wind farm with Snopt as the gradient-based optimizer.

# Arguments
- `final_layout_path::String`: path for the final layout file (.yaml)
- `opt_info_directory::String`: path to the directory where Snopt will store the .out files with the iteration history
- `layout_param::LayoutParamSet`: container holding all parameters related to the turbine starting layout
- `model_param::ModelParamSet`: container holding all parameters related to the turbine and farm analysis models
- `opt_algorithm::SnoptWECAlgorithm`: container holding all parameters for the Snopt+WEC optimizer
- `boundary::CircleBoundary`: container holding all parameters for the circular farm boundary
"""
function optimize_farm_layout(wind_farm_opt_with_TI, wind_farm_opt_no_TI, final_layout_path, opt_info_directory, layout_param, model_param, opt_algorithm::SnoptWECAlgorithm, boundary::CircleBoundary)


    #################################################################################
    # SET OPTIMIZATION PARAMETERS
    #################################################################################

    # get number of turbines in farm
    nturbines = length(layout_param.turbine_x)

    # set general lower and upper bounds for design variables
    lb = zeros(Int(nturbines*2)) .+ boundary.center[1] .- boundary.radius
    ub = zeros(Int(nturbines*2)) .+ boundary.center[2] .+ boundary.radius

    # get initial x value
    x = [deepcopy(layout_param.turbine_x);deepcopy(layout_param.turbine_y)]

    # initialize xopt array
    noptimizations = length(opt_algorithm.wec) + 1  # plus 1 b/c we are adding the starting point as the first column
    xopt_all = zeros(2*nturbines, noptimizations)
    xopt_all[:,1] = x

    # report initial objective value
    initial_aep, _, _, _, _ = wind_farm_opt_with_TI(x)
    println("starting objective value: ", -initial_aep*1e-6/opt_algorithm.objscale, " MWh\n")

    # set up options for SNOPT
    options = Dict{String, Any}()
    options["Derivative option"] = 1
    options["Major iteration limit"] = opt_algorithm.maxiter
    if opt_algorithm.checkgradients
        options["Verify level"] = 3
    else
        options["Verify level"] = -1
    end


    #################################################################################
    # RUN OPTIMIZATIONS
    #################################################################################

    # navigate to opt info directory (to store .out files from Snopt)
    original_directory = pwd()
    mkpath(opt_info_directory)
    cd(opt_info_directory)
    
    # run WEC optimizations
    for i = 1:length(opt_algorithm.wec)

        # set WEC value and strings for file paths
        if opt_algorithm.withTI[i]
            withTI_string = "TI-"
            model_param.farm_model_with_ti.wake_deficit_model.wec_factor[1] = opt_algorithm.wec[i]
            wec_string = string(round(model_param.farm_model_with_ti.wake_deficit_model.wec_factor[1],digits=1))
        else
            withTI_string = ""
            model_param.farm_model_no_ti.wake_deficit_model.wec_factor[1] = opt_algorithm.wec[i]
            wec_string = string(round(model_param.farm_model_no_ti.wake_deficit_model.wec_factor[1],digits=1))
        end
        println("Now running with WEC = ", wec_string, " and no local TI")

        # set Snopt options
        options["Summary file"] = "wec$(wec_string)-$(withTI_string)summary.out"
        options["Print file"] = "wec$(wec_string)-$(withTI_string)print.out"
        # options["Summary file"] = "summary.out"
        # options["Print file"] = "print.out"
        options["Major optimality tolerance"] = opt_algorithm.tol[i]

        # call Snopt
        if opt_algorithm.withTI[i]
            xopt_all[:,i+1], fopt, info = snopt(wind_farm_opt_with_TI, x, lb, ub, options) 
        else
            xopt_all[:,i+1], fopt, info = snopt(wind_farm_opt_no_TI, x, lb, ub, options)
        end
        println("Info: ", info, "\n")

    end

    # navigate back to original directory
    cd(original_directory)


    #################################################################################
    # DISPLAY/SAVE RESULTS
    #################################################################################

    # print optimization results
    final_aep, _, _, _, _ = wind_farm_opt_with_TI(xopt_all[:,end])
    println("end objective value: ", -final_aep*1e-6/opt_algorithm.objscale, " MWh\n")

    # # save WEC layout history to txt file
    # open("../results/final-layouts/" * main_file_path * "final-layout-wec-history-$(layout_number_string).txt", "w") do io
    #     write(io, "# initial, wec=3.0, wec=2.6, wec=2.2, wec=1.8, wec=1.4, wec=1.0, wec=1.0_withlocTI\n")
    #     writedlm(io, xopt_all)
    # end

    # extract final turbine locations
    turbine_x = copy(xopt_all[:,end][1:nturbines])
    turbine_y = copy(xopt_all[:,end][nturbines+1:end])

    # save final layout to YAML
    ff.write_turb_loc_YAML(final_layout_path, turbine_x, turbine_y; 
        title="Optimized turbine layout", titledescription="", 
        turbinefile="", locunits="m", wakemodelused="", windresourcefile="", aeptotal=final_aep*1e-6/opt_algorithm.objscale, 
        aepdirs="", aepunits="MWh", baseyaml="../data/initial-layouts/default_turbine_layout.yaml")

    return turbine_x, turbine_y
end



# function optimize_farm_layout(final_layout_path, opt_info_directory, layout_param, model_param, opt_algorithm::SnoptWECAlgorithm, boundary::CircleBoundary)


#     #################################################################################
#     # DEFINE OPTIMIZATION FUNCTIONS
#     #################################################################################

#     # set up objective wrapper function
#     function obj_with_TI(x)

#         # get number of turbines
#         nturbines = Int(length(x)/2)

#         # extract x and y locations of turbines from design variables vector
#         turbine_x = x[1:nturbines] 
#         turbine_y = x[nturbines+1:end]

#         # calculate AEP
#         AEP = opt_algorithm.objscale*ff.calculate_aep(turbine_x, turbine_y, layout_param.base_heights, model_param.turbine_design.rotor_diameter,
#                     model_param.turbine_design.hub_height, model_param.turbine_op.yaw, model_param.turbine_ct_models, model_param.turbine_design.generator_efficiency, model_param.turbine_design.cut_in_speed,
#                     model_param.turbine_design.cut_out_speed, model_param.turbine_design.rated_speed, model_param.turbine_design.rated_power, model_param.wind_resource_model, model_param.turbine_power_models, model_param.farm_model_with_ti,
#                     rotor_sample_points_y=model_param.velocity_sampling.rotor_sample_points_y, rotor_sample_points_z=model_param.velocity_sampling.rotor_sample_points_z)
        
#         # return the objective as an array
#         return [AEP]
#     end

#     function obj_no_TI(x)

#         # get number of turbines
#         nturbines = Int(length(x)/2)

#         # extract x and y locations of turbines from design variables vector
#         turbine_x = x[1:nturbines] 
#         turbine_y = x[nturbines+1:end]

#         # calculate AEP
#         AEP = opt_algorithm.objscale*ff.calculate_aep(turbine_x, turbine_y, layout_param.base_heights, model_param.turbine_design.rotor_diameter,
#                     model_param.turbine_design.hub_height, model_param.turbine_op.yaw, model_param.turbine_ct_models, model_param.turbine_design.generator_efficiency, model_param.turbine_design.cut_in_speed,
#                     model_param.turbine_design.cut_out_speed, model_param.turbine_design.rated_speed, model_param.turbine_design.rated_power, model_param.wind_resource_model, model_param.turbine_power_models, model_param.farm_model_no_ti,
#                     rotor_sample_points_y=model_param.velocity_sampling.rotor_sample_points_y, rotor_sample_points_z=model_param.velocity_sampling.rotor_sample_points_z)

#         # return the objective as an array
#         return [AEP]
#     end

#     # set up constraint wrapper function
#     function con(x)

#         # get number of turbines
#         nturbines = Int(length(x)/2)
        
#         # extract x and y locations of turbines from design variables vector
#         turbine_x = x[1:nturbines]
#         turbine_y = x[nturbines+1:end]

#         # get constraint values
#         spacing_con = 2.0*model_param.turbine_design.rotor_diameter[1] .- ff.turbine_spacing(turbine_x,turbine_y)
#         boundary_con = ff.circle_boundary(boundary.center, boundary.radius, turbine_x, turbine_y)
        
#         return [spacing_con; boundary_con]
#     end

#     # set up optimization problem wrapper function
#     function wind_farm_opt_with_TI(x)

#         # calculate the objective and jacobian (negative sign in order to maximize AEP)
#         AEP = -obj_with_TI(x)[1]
#         dAEP_dx = -ForwardDiff.jacobian(obj_with_TI, x)

#         # calculate constraint and jacobian
#         c = con(x)
#         dc_dx = ForwardDiff.jacobian(con, x)

#         # set fail flag to false
#         fail = false

#         # return objective, constraint, and jacobian values
#         return AEP, c, dAEP_dx, dc_dx, fail
#     end

#     function wind_farm_opt_no_TI(x)

#         # calculate the objective and jacobian (negative sign in order to maximize AEP)
#         AEP = -obj_no_TI(x)[1]
#         dAEP_dx = -ForwardDiff.jacobian(obj_no_TI, x)

#         # calculate constraint and jacobian
#         c = con(x)
#         dc_dx = ForwardDiff.jacobian(con, x)

#         # set fail flag to false
#         fail = false

#         # return objective, constraint, and jacobian values
#         return AEP, c, dAEP_dx, dc_dx, fail
#     end


#     #################################################################################
#     # SET OPTIMIZATION PARAMETERS
#     #################################################################################

#     # get number of turbines in farm
#     nturbines = length(layout_param.turbine_x)

#     # set general lower and upper bounds for design variables
#     lb = zeros(Int(nturbines*2)) .+ boundary.center[1] .- boundary.radius
#     ub = zeros(Int(nturbines*2)) .+ boundary.center[2] .+ boundary.radius

#     # get initial x value
#     x = [deepcopy(layout_param.turbine_x);deepcopy(layout_param.turbine_y)]

#     # initialize xopt array
#     noptimizations = length(opt_algorithm.wec) + 1  # plus 1 b/c we are adding the starting point as the first column
#     xopt_all = zeros(2*nturbines, noptimizations)
#     xopt_all[:,1] = x

#     # report initial objective value
#     println("starting objective value: ", obj_with_TI(x)[1], "\n")

#     # set up options for SNOPT
#     options = Dict{String, Any}()
#     options["Derivative option"] = 1
#     options["Major iteration limit"] = opt_algorithm.maxiter
#     if opt_algorithm.checkgradients
#         options["Verify level"] = 3
#     else
#         options["Verify level"] = -1
#     end


#     #################################################################################
#     # RUN OPTIMIZATIONS
#     #################################################################################

#     # navigate to opt info directory (to store .out files from Snopt)
#     original_directory = pwd()
#     mkpath(opt_info_directory)
#     cd(opt_info_directory)
    
#     # run WEC optimizations
#     for i = 1:length(opt_algorithm.wec) 

#         # set WEC value and strings for file paths
#         if opt_algorithm.withTI[i]
#             withTI_string = "TI-"
#             model_param.farm_model_with_ti.wake_deficit_model.wec_factor[1] = opt_algorithm.wec[i]
#             wec_string = string(round(model_param.farm_model_with_ti.wake_deficit_model.wec_factor[1],digits=1))
#         else
#             withTI_string = ""
#             model_param.farm_model_no_ti.wake_deficit_model.wec_factor[1] = opt_algorithm.wec[i]
#             wec_string = string(round(model_param.farm_model_no_ti.wake_deficit_model.wec_factor[1],digits=1))
#         end
#         println("Now running with WEC = ", wec_string, " and no local TI")

#         # set Snopt options
#         options["Summary file"] = "wec$(wec_string)-$(withTI_string)summary.out"
#         options["Print file"] = "wec$(wec_string)-$(withTI_string)print.out"
#         options["Major optimality tolerance"] = opt_algorithm.tol[i]

#         # call Snopt
#         if opt_algorithm.withTI[i]
#             xopt_all[:,i+1], fopt, info = snopt(wind_farm_opt_with_TI, x, lb, ub, options) 
#         else
#             xopt_all[:,i+1], fopt, info = snopt(wind_farm_opt_no_TI, x, lb, ub, options)
#         end
#         println("Info: ", info, "\n")

#     end

#     # navigate back to original directory
#     cd(original_directory)


#     #################################################################################
#     # DISPLAY/SAVE RESULTS
#     #################################################################################

#     # print optimization results
#     final_obj = obj_with_TI(xopt_all[:,end])[1]
#     println("end objective value: ", final_obj)

#     # # save WEC layout history to txt file
#     # open("../results/final-layouts/" * main_file_path * "final-layout-wec-history-$(layout_number_string).txt", "w") do io
#     #     write(io, "# initial, wec=3.0, wec=2.6, wec=2.2, wec=1.8, wec=1.4, wec=1.0, wec=1.0_withlocTI\n")
#     #     writedlm(io, xopt_all)
#     # end

#     # extract final turbine locations
#     turbine_x = copy(xopt_all[:,end][1:nturbines])
#     turbine_y = copy(xopt_all[:,end][nturbines+1:end])

#     # save final layout to YAML
#     ff.write_turb_loc_YAML(final_layout_path, turbine_x, turbine_y; 
#         title="Optimized turbine layout", titledescription="", 
#         turbinefile="", locunits="m", wakemodelused="", windresourcefile="", aeptotal=final_obj*1e-6/opt_algorithm.objscale, 
#         aepdirs="", aepunits="MWh", baseyaml="../data/initial-layouts/default_turbine_layout.yaml")

#     return turbine_x, turbine_y
# end



"""
    optimize_farm_layout(final_layout_path, opt_info_directory, layout_param, model_param, opt_algorithm::SnoptWECAlgorithm, boundary::FreeFormBoundary)

Runs the turbine layout optimization for a circular boundary wind farm with Snopt as the gradient-based optimizer.

# Arguments
- `final_layout_path::String`: path for the final layout file (.yaml)
- `opt_info_directory::String`: path to the directory where Snopt will store the .out files with the iteration history
- `layout_param::LayoutParamSet`: container holding all parameters related to the turbine starting layout
- `model_param::ModelParamSet`: container holding all parameters related to the turbine and farm analysis models
- `opt_algorithm::SnoptWECAlgorithm`: container holding all parameters for the Snopt+WEC optimizer
- `boundary::CircleBoundary`: container holding all parameters for the circular farm boundary
"""
function optimize_farm_layout(wind_farm_opt_with_TI, wind_farm_opt_no_TI, final_layout_path, opt_info_directory, layout_param, model_param, opt_algorithm::SnoptWECAlgorithm, boundary::FreeFormBoundary)


    #################################################################################
    # SET OPTIMIZATION PARAMETERS
    #################################################################################

    # get number of turbines in farm
    nturbines = length(layout_param.turbine_x)

    # set general lower and upper bounds for design variables
    lb = zeros(Int(nturbines*2)) .+ minimum(boundary.vertices)
    ub = zeros(Int(nturbines*2)) .+ maximum(boundary.vertices)

    # get initial x value
    x = [deepcopy(layout_param.turbine_x);deepcopy(layout_param.turbine_y)]

    # initialize xopt array
    noptimizations = length(opt_algorithm.wec) + 1  # plus 1 b/c we are adding the starting point as the first column
    xopt_all = zeros(2*nturbines, noptimizations)
    xopt_all[:,1] = x

    # report initial objective value
    initial_aep, _, _, _, _ = wind_farm_opt_with_TI(x)
    println("starting objective value: ", -initial_aep*1e-6/opt_algorithm.objscale, " MWh\n")

    # set up options for SNOPT
    options = Dict{String, Any}()
    options["Derivative option"] = 1
    options["Major iteration limit"] = opt_algorithm.maxiter
    if opt_algorithm.checkgradients
        options["Verify level"] = 3
    else
        options["Verify level"] = -1
    end


    #################################################################################
    # RUN OPTIMIZATIONS
    #################################################################################

    # navigate to opt info directory (to store .out files from Snopt)
    original_directory = pwd()
    mkpath(opt_info_directory)
    cd(opt_info_directory)
    
    # run WEC optimizations
    for i = 1:length(opt_algorithm.wec)

        # set WEC value and strings for file paths
        if opt_algorithm.withTI[i]
            withTI_string = "TI-"
            model_param.farm_model_with_ti.wake_deficit_model.wec_factor[1] = opt_algorithm.wec[i]
            wec_string = string(round(model_param.farm_model_with_ti.wake_deficit_model.wec_factor[1],digits=1))
        else
            withTI_string = ""
            model_param.farm_model_no_ti.wake_deficit_model.wec_factor[1] = opt_algorithm.wec[i]
            wec_string = string(round(model_param.farm_model_no_ti.wake_deficit_model.wec_factor[1],digits=1))
        end
        println("Now running with WEC = ", wec_string, " and no local TI")

        # set Snopt options
        options["Summary file"] = "wec$(wec_string)-$(withTI_string)summary.out"
        options["Print file"] = "wec$(wec_string)-$(withTI_string)print.out"
        options["Major optimality tolerance"] = opt_algorithm.tol[i]

        # call Snopt
        if opt_algorithm.withTI[i]
            xopt_all[:,i+1], fopt, info = snopt(wind_farm_opt_with_TI, x, lb, ub, options) 
        else
            xopt_all[:,i+1], fopt, info = snopt(wind_farm_opt_no_TI, x, lb, ub, options)
        end
        println("Info: ", info, "\n")

    end

    # navigate back to original directory
    cd(original_directory)


    #################################################################################
    # DISPLAY/SAVE RESULTS
    #################################################################################

    # print optimization results
    final_aep, _, _, _, _ = wind_farm_opt_with_TI(xopt_all[:,end])
    println("end objective value: ", -final_aep*1e-6/opt_algorithm.objscale, " MWh\n")

    # # save WEC layout history to txt file
    # open("../results/final-layouts/" * main_file_path * "final-layout-wec-history-$(layout_number_string).txt", "w") do io
    #     write(io, "# initial, wec=3.0, wec=2.6, wec=2.2, wec=1.8, wec=1.4, wec=1.0, wec=1.0_withlocTI\n")
    #     writedlm(io, xopt_all)
    # end

    # extract final turbine locations
    turbine_x = copy(xopt_all[:,end][1:nturbines])
    turbine_y = copy(xopt_all[:,end][nturbines+1:end])

    # save final layout to YAML
    ff.write_turb_loc_YAML(final_layout_path, turbine_x, turbine_y; 
        title="Optimized turbine layout", titledescription="", 
        turbinefile="", locunits="m", wakemodelused="", windresourcefile="", aeptotal=final_aep*1e-6/opt_algorithm.objscale, 
        aepdirs="", aepunits="MWh", baseyaml="../data/initial-layouts/default_turbine_layout.yaml")

    return turbine_x, turbine_y
end
