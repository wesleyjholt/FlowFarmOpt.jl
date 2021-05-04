abstract type WindRose end 


"""
    NantucketWindRose(ndirs)

Container holding parameters for the Nantucket wind rose

# Arguments
- `ndirs::Float64`: number of wind direction bins
"""
struct NantucketWindRose{I, F} <: WindRose
    ndirs::I
    nspeeds::I
    air_density::F
    ambient_ti::F
    shear_exponent::F
end
NantucketWindRose() = NantucketWindRose(10, 1, 1.225, 0.108, 0.31)
NantucketWindRose(ndirs) = NantucketWindRose(ndirs,  1, 1.225, 0.108, 0.31)
NantucketWindRose(ndirs, nspeeds) = NantucketWindRose(ndirs, nspeeds, 1.225, 0.108, 0.31)

"""
    HornsRevWindRose(ndirs, nspeeds)

Container holding parameters for the Horns Rev wind rose

# Arguments
- `ndirs::Float64`: number of wind direction bins
- `nspeeds::Float64`: number of wind speed bins
"""
struct HornsRevWindRose{I, F} <: WindRose
    ndirs::I
    nspeeds::I
    air_density::F
    ambient_ti::F
    shear_exponent::F
end
HornsRevWindRose() = HornsRevWindRose(10, 1, 1.225, 0.108, 0.31)
HornsRevWindRose(ndirs) = HornsRevWindRose(ndirs, 1, 1.225, 0.108, 0.31)
HornsRevWindRose(ndirs, nspeeds) = HornsRevWindRose(ndirs, nspeeds, 1.225, 0.108, 0.31)


"""
    resample_wind_resource(windrose::NantucketWindRose)

Resamples the Nantucket wind resource for the specified number of directions. (average speeds)

# Arguments
-`windrose::NantucketWindRose`: container holding parameters for the Nantucket wind rose
"""
function resample_wind_resource(windrose::NantucketWindRose)

    # load original wind resource
    input_windrose_file_path = "../data/windrose-files/nantucket/nantucket-windrose-ave-speeds.txt"
    windrose_data = readdlm(input_windrose_file_path, skipstart=1)
    directions_orig = windrose_data[:,1]
    speeds_orig = windrose_data[:,2]
    direction_probabilities_orig = windrose_data[:,3]

    # resample discretized wind resource
    directions, direction_probabilities, speeds = resample_discretized_windrose_average_speeds(directions_orig, direction_probabilities_orig, speeds_orig, windrose.ndirs)

    # write to a txt file
    output_windrose_file_path = "../data/windrose-files/nantucket/nantucket-windrose-ave-speeds-$(lpad(windrose.ndirs, 3, "0"))dirs.txt"
    open(output_windrose_file_path, "w") do io
        write(io, "# directions, average speeds, frequencies\n")
        writedlm(io, [directions speeds direction_probabilities])
    end
end

"""
    resample_wind_resource(windrose::HornsRevWindRose)

Resamples the Horns Rev 1 wind resource for the specified number of directions and speeds.
If only one speed is asked for, then a wind rose with average speeds for direction is created.

# Arguments
- `windrose::HornsRevWindRose`: container holding parameters for the Horns Rev 1 wind rose
"""
function resample_wind_resource(windrose::HornsRevWindRose)
    # set original values for direction and speed joint distribution (from paper by Ju Feng and Wen Zhong Shen, 2015, "Modelling Wind for Wind Farm Layout Optimization Using Joint Distribution of Wind Speed and Wind Direction")
    directions_orig = [0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0, 210.0, 240.0, 270.0, 300.0, 330.0]
    direction_probabilities_orig = [0.0482, 0.0406, 0.0359, 0.0527, 0.0912, 0.0697, 0.0917, 0.1184, 0.1241, 0.1134, 0.117, 0.0969]
    speed_weibull_params_orig = [2.09 8.89; 2.13 9.27; 2.29 8.23; 2.30 9.78; 2.67 11.64; 2.45 11.03; 2.51 11.50; 2.40 11.92; 2.35 11.49; 2.27 11.08; 2.24 11.34; 2.19 10.76]

    # resample discretized wind rose
    if windrose.nspeeds != 1
        # wind speed discretized distributions
        directions, direction_probabilities, speeds, speed_probabilities = resample_discretized_windrose_multiple_speeds_weibull(directions_orig, direction_probabilities_orig, speed_weibull_params_orig, windrose.ndirs, windrose.nspeeds)
    else
        # average wind speeds
        directions, direction_probabilities, speeds = resample_discretized_windrose_average_speeds_weibull(directions_orig, direction_probabilities_orig, speed_weibull_params_orig, windrose.ndirs)
        speed_probabilities = []
    end

    # write to a yaml file
    output_windrose_filepath = "../data/windrose-files/horns-rev/horns-rev-windrose-$(lpad(windrose.ndirs,3,"0"))dirs-$(lpad(windrose.nspeeds,2,"0"))speeds.yaml"
    title = "Horns Rev 1 Wind Rose, $(windrose.ndirs) wind directions, $(windrose.nspeeds) wind speeds"
    description = "Wind resource conditions, using direction and speed distribution from the paper by Ju Feng and Wen Zhong Shen: \nModelling Wind for Wind Farm Layout Optimization Using Joint Distribution of Wind Speed and Wind Direction \nhttps://backend.orbit.dtu.dk/ws/portalfiles/portal/110639016/Modelling_Wind_for_Wind_Farm_Layout.pdf"
    turbulence_intensity = 0.75
    write_windrose_yaml(output_windrose_filepath, directions, direction_probabilities, speeds, speed_probabilities, turbulence_intensity, title=title, description=description)
end


"""
    get_wind_resource_model(windrose::NantucketWindRose)

Creates a FlowFarm wind resource model for the Nantucket wind resource.

# Arguments
- `windrose::NantucketWindRose`: container holding parameters for the Nantucket wind rose
"""
function get_wind_resource_model(windrose::NantucketWindRose)
    # get wind directions and speeds joint pmf
    windrose_file_path = "../data/windrose-files/nantucket/nantucket-windrose-ave-speeds-$(lpad(windrose.ndirs, 3, "0"))dirs.txt"
    windrose_data = readdlm(windrose_file_path, skipstart=1)
    wind_directions = windrose_data[:,1]*pi/180          # radians
    wind_speeds = windrose_data[:,2]                     # m/s
    wind_probabilities = windrose_data[:,3]
    
    # set remaining wind parameters
    nstates = length(wind_speeds)
    ambient_tis = zeros(nstates) .+ windrose.ambient_ti
    measurement_height = zeros(nstates) .+ 80.0

    # set FlowFarm wind resource models
    wind_shear_model = ff.PowerLawWindShear(windrose.shear_exponent)
    wind_resource = ff.DiscretizedWindResource(wind_directions, wind_speeds, wind_probabilities, measurement_height, windrose.air_density, ambient_tis, wind_shear_model)
    
    return wind_resource
end

"""
    get_wind_resource_model(windrose::HornsRevWindRose)

Creates a FlowFarm wind resource model for the Horns Rev 1 wind resource.

# Arguments
- `windrose::HornsRevWindRose`: container holding parameters for the Horns Rev 1 wind rose
"""
function get_wind_resource_model(windrose::HornsRevWindRose)

    # get wind directions and speeds joint pmf
    windrose_file_path = "../data/windrose-files/horns-rev/horns-rev-windrose-$(lpad(windrose.ndirs,3,"0"))dirs-$(lpad(windrose.nspeeds,2,"0"))speeds.yaml"
    if windrose.nspeeds != 1
        # wind speed discretized distributions
        wind_directions, wind_speeds, wind_probabilities, ambient_ti = ff.get_wind_rose_YAML(windrose_file_path)
    else
        # average wind speeds
        wind_data = YAML.load_file(windrose_file_path)["definitions"]["wind_inflow"]["properties"]
        wind_directions = wind_data["direction"]["bins"]
        wind_speeds = wind_data["speed"]["bins"]
        wind_probabilities = wind_data["direction"]["frequency"]
        ambient_ti = wind_data["turbulence_intensity"]["default"]
    end
    wind_directions *= pi/180.0
    
    # set remaining wind parameters
    nstates = length(wind_speeds)
    ambient_tis = zeros(nstates) .+ windrose.ambient_ti
    measurement_height = zeros(nstates) .+ 80.0

    # set FlowFarm wind resource models
    wind_shear_model = ff.PowerLawWindShear(windrose.shear_exponent)
    wind_resource = ff.DiscretizedWindResource(wind_directions, wind_speeds, wind_probabilities, measurement_height, windrose.air_density, ambient_tis, wind_shear_model)
    
    return wind_resource
end



function resample_discretized_windrose_average_speeds(directions_orig, direction_probabilities_orig, speeds_orig, ndirs)

    # create cyclical data
    directions_orig_cyclical = [directions_orig[end-1:end] .- 360.0; directions_orig; directions_orig[1:2] .+ 360.0]
    direction_probabilities_orig_cyclical = [direction_probabilities_orig[end-1:end]; direction_probabilities_orig; direction_probabilities_orig[1:2]]
    speeds_orig_cyclical = [speeds_orig[end-1:end]; speeds_orig; speeds_orig[1:2]]

    # get direction vector
    directions = range(0, stop=360*(ndirs-1)/ndirs, length=ndirs)

    # use akima spline to interpolate the pmf to the specified number of wind directions
    direction_probabilities = akima(directions_orig_cyclical, direction_probabilities_orig_cyclical, directions)
    direction_probabilities /= sum(direction_probabilities)

    # use akima spline to interpolate speeds
    speeds = akima(directions_orig_cyclical, speeds_orig_cyclical, directions)

    return directions, direction_probabilities, speeds
end



function resample_discretized_windrose_average_speeds_weibull(directions_orig, direction_probabilities_orig, speed_weibull_params_orig, ndirs)

    # create cyclical data
    directions_orig_cyclical = [directions_orig[end-1:end] .- 360.0; directions_orig; directions_orig[1:2] .+ 360.0]
    direction_probabilities_orig_cyclical = [direction_probabilities_orig[end-1:end]; direction_probabilities_orig; direction_probabilities_orig[1:2]]
    speed_weibull_params_orig_cyclical = [speed_weibull_params_orig[end-1:end,:]; speed_weibull_params_orig; speed_weibull_params_orig[1:2,:]]

    # get direction vector
    directions = range(0, stop=360*(ndirs-1)/ndirs, length=ndirs)

    # use akima spline to interpolate the direction probabilities and normalize the interpolated probabilities to make them an actual probability mass function 
    direction_probabilities = akima(directions_orig_cyclical, direction_probabilities_orig_cyclical, directions)
    direction_probabilities /= sum(direction_probabilities)

    # use akima spline to interpolate the speed weibull parameters to the specified number of wind directions
    speed_weibull_params = zeros(ndirs, 2)
    speed_weibull_params[:,1] = akima(directions_orig_cyclical, speed_weibull_params_orig_cyclical[:,1], directions)
    speed_weibull_params[:,2] = akima(directions_orig_cyclical, speed_weibull_params_orig_cyclical[:,2], directions)

    # get average speeds
    average_speeds = zeros(ndirs)
    n_integration_bins = 200
    delta_speed_integration = 25.0/n_integration_bins
    speed_integration_bins = range(0.0 + delta_speed_integration/2, stop = 25.0 - delta_speed_integration/2, step = delta_speed_integration)
    for i = 1:ndirs
        dist = Weibull(speed_weibull_params[i,:]...)
        average_speeds[i] = sum(delta_speed_integration .* speed_integration_bins .* pdf.(dist, speed_integration_bins))
    end

    return directions, direction_probabilities, average_speeds
end



function resample_discretized_windrose_multiple_speeds_weibull(directions_orig, direction_probabilities_orig, speed_weibull_params_orig, ndirs, nspeeds)

    # create cyclical data
    directions_orig_cyclical = [directions_orig[end-1:end] .- 360.0; directions_orig; directions_orig[1:2] .+ 360.0]
    direction_probabilities_orig_cyclical = [direction_probabilities_orig[end-1:end]; direction_probabilities_orig; direction_probabilities_orig[1:2]]
    speed_weibull_params_orig_cyclical = [speed_weibull_params_orig[end-1:end,:]; speed_weibull_params_orig; speed_weibull_params_orig[1:2,:]]

    # get direction and speed vectors
    directions = range(0, stop=360*(ndirs-1)/ndirs, length=ndirs)
    delta_speed = 25.0/nspeeds
    speeds = range(0.0 + delta_speed/2, stop = 25.0 - delta_speed/2, step = delta_speed)
    speed_bin_edges = range(0.0, stop=25.0, length=nspeeds+1)

    # use akima spline to interpolate the probabilities and weibull parameters to the specified number of wind directions
    direction_probabilities = akima(directions_orig_cyclical, direction_probabilities_orig_cyclical, directions)
    speed_weibull_params = zeros(ndirs, 2)
    speed_weibull_params[:,1] = akima(directions_orig_cyclical, speed_weibull_params_orig_cyclical[:,1], directions)
    speed_weibull_params[:,2] = akima(directions_orig_cyclical, speed_weibull_params_orig_cyclical[:,2], directions)

    # normalize the interpolated probabilities to make them an actual probability mass function 
    direction_probabilities /= sum(direction_probabilities)

    # get the speed probabilities for each direction
    speed_probabilities = fill(zeros(length(speeds)), ndirs)
    for i = 1:ndirs
        dist = Weibull(speed_weibull_params[i,:]...)
        speed_probabilities[i] = cdf.(dist, speed_bin_edges[2:end]) - cdf.(dist, speed_bin_edges[1:end-1])
        # normalize
        speed_probabilities[i] ./= sum(speed_probabilities[i])
    end

    return directions, direction_probabilities, speeds, speed_probabilities
end


"""
    write_windrose_yaml(output_windrose_filepath, directions, direction_probabilities, speeds, speed_probabilities, turbulence_intensity; title="", description="", template_file="../data/windrose-files/windrose-template.yaml")

Saves wind resource information to a YAML file.

# Arguments
- `output_windrose_filepath::String`: file path for output wind resource
- `directions::Array{Float64,1}`: contains wind direction bins
- `direction_probabilities::Array{Float64,1}`: contains wind direction probabilities
- `speeds::Array{Float64,1}`: contains wind speed bins
- `speed_probabilities::Array{Array{Float64,1},1}`: contains wind speed probabilities (one array for each direction)
- `turbulence_intensity::Float64`: turbulence intensity
"""
function write_windrose_yaml(output_windrose_filepath, directions, direction_probabilities, speeds, speed_probabilities, turbulence_intensity; title="", description="", template_file = "../data/windrose-files/windrose-template.yaml")

    # save wind data to a data structure (from a template)
    data = YAML.load_file(template_file)
    data["title"] = title
    data["description"] = description
    wind_data = data["definitions"]["wind_inflow"]["properties"]
    wind_data["speed"]["bins"] = speeds
    wind_data["speed"]["frequency"] = speed_probabilities
    wind_data["speed"]["minimum"] = minimum(speeds)
    wind_data["speed"]["maximum"] = maximum(speeds)
    wind_data["direction"]["bins"] = directions
    wind_data["direction"]["frequency"] = direction_probabilities
    wind_data["turbulence_intensity"]["default"] = turbulence_intensity

    # write to a yaml file
    write = YAML.write_file(output_windrose_filepath, data)

end