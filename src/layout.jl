"""
    LayoutParamSet(turbine_x, turbine_y, base_heights)

Container holding coordinates and base heights for the turbines in a farm.

# Arguments
- `turbine_x::Array{Float64,1}`: x-coordinates for the turbines in the farm
- `turbine_y::Array{Float64,1}`: y-coordinates for the turbines in the farm
- `base_heights::Array{Float64,1}`: base heights for the turbines (not including the tower heights)
"""
mutable struct LayoutParamSet{}
    turbine_x
    turbine_y
    base_heights
    aep
end

"""
    layout_set(initial_layout_path)

Creates a LayoutParamSet object to hold the turbine layout for a farm. Assumes the entire farm is on level ground.

# Arguments
- `initial_layout_path::String`: path to the text file with the turbine coordinates
"""
function layout_set(layout_path; file_type="text")
    if file_type=="text"
        turbine_locations = readdlm(layout_path, skipstart=1)
        turbine_x = turbine_locations[:,1]
        turbine_y = turbine_locations[:,2]
        aep = ""
    elseif file_type=="yaml"
        layout_data = YAML.load_file(layout_path)
        nturbines = length(layout_data["definitions"]["position"]["items"])
        turbine_x = zeros(nturbines)
        turbine_y = zeros(nturbines)
        for j = 1:nturbines
            turbine_x[j] = layout_data["definitions"]["position"]["items"][j][1]
            turbine_y[j] = layout_data["definitions"]["position"]["items"][j][2]
        end
        aep = layout_data["definitions"]["plant_energy"]["properties"]["annual_energy_production"]["default"]
    else
        error("Invalid file type.")
    end
    base_heights = zeros(length(turbine_x))
    layout_param = LayoutParamSet(turbine_x, turbine_y, base_heights, aep)
    return layout_param
end


