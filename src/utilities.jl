# using Snopt

function filter_aeps(aeps_raw; layouts_must_match=true)
    n_ndirs = length(aeps_raw[:,1])
    n_nspeeds = length(aeps_raw[1,:])
    if layouts_must_match
        n_layouts = maximum(length.(aeps_raw))
        good_indices = trues(n_layouts)
        for k = 1:n_layouts
            for i = 1:n_ndirs
                for j = 1:n_nspeeds
                    if isnan(aeps_raw[i,j][k])
                        good_indices[k] = false
                    end
                end
            end
        end
        aeps_filtered = fill(zeros(0), n_ndirs, n_nspeeds)
        for i = 1:n_ndirs
            for j = 1:n_nspeeds
                aeps_filtered[i,j] = aeps_raw[i,j][good_indices]
            end
        end
    else
        n_layouts = maximum(length.(aeps_raw))
        aeps_filtered = deepcopy(aeps_raw)
        for k = 1:n_layouts
            for i = 1:n_ndirs
                for j = 1:n_nspeeds
                    if length(aeps_filtered[i,j]) >= k
                        if isnan(aeps_filtered[i,j][k])
                            splice!(aeps_filtered[i,j],k)
                        end
                    end
                end
            end
        end
    end

    return aeps_filtered
end


function normalize_aeps(aeps::Array{Array{Float64,1},2})
    n_ndirs = length(aeps[:,1])
    n_nspeeds = length(aeps[1,:])
    n_layouts = minimum(length.(aeps))
    aeps_normalized = fill(zeros(0), n_ndirs, n_nspeeds)
    for i = 1:n_ndirs
        for j = 1:n_nspeeds
            aeps_normalized[i,j] = aeps[i,j] ./ aeps[end,end]
        end
    end
    return aeps_normalized
end


function get_aep_values_from_file_names(data_file_name::String)
    # get the AEP values from the files
    aeps = readdlm(data_file_name, skipstart=1)[:,2]
    return aeps
end


function get_aep_values_from_file_names(data_file_names::Array{String,1})
    # get the AEP values from the files
    ngroups = length(data_file_names)
    aeps = fill(zeros(1), ngroups)
    for i = 1:ngroups
        aeps[i] = readdlm(data_file_names[i], skipstart=1)[:,2]
    end
    return aeps
end


function get_aep_values_from_file_names(data_file_names::Array{String,2})
    # specify the number of groups
    n_ndirs_vec = length(data_file_names[:,1])
    n_nspeeds_vec = length(data_file_names[1,:])
    # get AEP values from files
    aeps = fill(zeros(1), n_ndirs_vec, n_nspeeds_vec)
    for i = 1:n_ndirs_vec
        for j = 1:n_nspeeds_vec
            aeps[i,j] = readdlm(data_file_names[i,j], skipstart=1)[:,2]
        end
    end
    return aeps
end


function get_log_spaced_montecarlo_points(n_points, n_samples)

    function obj(x)
        -x[1]/1e6
    end

    function con(x, n)
        x_vec = get_x(x[1])
        c = sum(10 .^ x_vec)/n - 1
    end

    function opt(x, n)

        f = obj(x)
        # dfdx = ForwardDiff.jacobian(obj, x)

        c = con(x, n)
        # dcdx = ForwardDiff.jacobian(con, x)

        fail = false

        return f, c, fail
    end

    function get_x(x_scalar, n_points)
        x_vec = zeros(n_points)
        for i = 1:n_points
            x_vec[i] = 1.0 + x_scalar*(i-1)
        end
        return x_vec
    end

    x = [0.5]
    lb = zeros(n_points)
    ub = zeros(n_points) .+ 2.0
    options = Dict{String, Any}()
    options["Derivative option"] = 0

    opt(x) = opt(x, n_samples)
    get_x(x) = get_x(x, n_points)

    xopt, fopt, info = snopt(opt, x, lb, ub, options)

    xopt_vec = get_x(xopt[1], n_points)

    return Int.(floor.(10 .^ xopt_vec))

end


function calc_mean_upper10percent(x)
    n_x = length(x)
    index_90percentile = Int(round(0.8*n_x))
    sorted_indices = sortperm(x)
    mean_upper10percent_aeps = mean(x[sorted_indices[index_90percentile:end]])
    return mean_upper10percent_aeps
end

function calc_std_upper10percent(x)
    n_x = length(x)
    index_90percentile = Int(round(0.8*n_x))
    sorted_indices = sortperm(x)
    std_upper10percent_aeps = std(x[sorted_indices[index_90percentile:end]])
    return std_upper10percent_aeps
end
