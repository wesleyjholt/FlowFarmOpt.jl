"""
    FarmParams(turbine_x, turbine_y, base_heights, turbine_types)

Container to hold the initial random layout of the turbines

# Arguments
- `turbine_x::Array{Float64,1}`: turbine x coordinates
- `turbine_y::Array{Float64,1}`: turbine y coordinates
- `base_heights::Array{Float64,1}`: turbine base heights
- `turbine_type::Array{Int64,1}`: turbine types
"""
struct Layout{AF, S, I}
    turbine_x::AF
    turbine_y::AF
    layout_number::S
    nturbines::I
    base_heights::AF
end
Layout(x, y) = Layout(x, y, "000", length(x), zeros(length(x)))
Layout(x, y, a) = Layout(x, y, a, length(x), zeros(length(x)))


"""
    VelocitySampling(rotor_points_y, rotor_points_z)

Container holding rotor sampling points

# Arguments
- `rotor_sample_points_y::Array{Float64,1}`
- `rotor_sample_points_y::Array{Float64,1}`
"""
struct RotorSampling{AF}
    rotor_sample_points_y::AF
    rotor_sample_points_z::AF
end
RotorSampling() = RotorSampling([0.0], [0.0])
