abstract type FarmBoundary end

"""
    CircleBoundary(radius, center)

Container for parameters describing a circular wind farm boundary.

# Arguments
- `radius::Float64`: radius of the circular farm boundary
- `tol::Array{Float64,1}`: center of the farm
"""
struct CircleBoundary{F, AF} <: FarmBoundary
    radius::F
    center::AF
end
CircleBoundary() = CircleBoundary(1000.0, [0.0, 0.0])
CircleBoundary(a) = CircleBoundary(Float64(a), [0.0, 0.0])


"""
    RectangleBoundary(base, height, center)

Container for parameters describing a rectangular wind farm boundary.

# Arguments
- `base::Float64`: base of the rectangular farm boundary
- `height::Float64`: height of the rectangular farm boundary
- `tol::Array{Float64,1}`: center of the farm
"""
struct RectangleBoundary{F, AF, AF2} <: FarmBoundary
    base::F
    height::F
    center::AF
    vertices::AF2
    normals::AF2
end
RectangleBoundary() = RectangleBoundary(1000.0, 1000.0, [0.0, 0.0], get_rectangle_boundary_vertices(1000.0, 1000.0, [0.0, 0.0]), get_rectangle_boundary_normals(1000.0, 1000.0, [0.0, 0.0]))
RectangleBoundary(s) = RectangleBoundary(Float64(s), Float64(s), [0.0, 0.0], get_rectangle_boundary_vertices(s, s, [0.0, 0.0]), get_rectangle_boundary_normals(s, s, [0.0, 0.0]))
RectangleBoundary(b, h) = RectangleBoundary(Float64(b), Float64(h), [0.0, 0.0], get_rectangle_boundary_vertices(b, h, [0.0, 0.0]), get_rectangle_boundary_normals(b, h, [0.0, 0.0]))


# """
#     ParallelogramBoundary(base, height, angle, center)

# Container for parameters describing a parallelogram-shaped wind farm boundary.

# # Arguments
# - `base::Float64`: base of the parallelogram-shaped farm boundary
# - `height::Float64`: height of the parallelogram-shaped farm boundary
# - `angle::Float64`: angle of stretch (in degrees)
# - `tol::Array{Float64,1}`: center of the farm
# """
# struct ParallelogramBoundary{F, AF} <: FarmBoundary
#     base::F
#     height::F
#     angle::F
#     center::AF
# end
# ParallelogramBoundary() = ParallelogramBoundary(1000.0, 1000.0, 10.0, [0.0, 0.0])
# ParallelogramBoundary(a, b) = ParallelogramBoundary(a, b, 10.0, [0.0, 0.0])
# ParallelogramBoundary(a, b, c) = ParallelogramBoundary(a, b, c, [0.0, 0.0])


"""
    FreeFormBoundary(vertices)

Container for parameters describing a free-form wind farm boundary.

# Arguments
- `vertices::Array{Float64,2}`: coordinates (x,y) for the vertices of the free-form boundary (counter-clockwise)
"""
struct FreeFormBoundary{AF} <: FarmBoundary
    vertices::AF
    normals::AF
end
function FreeFormBoundary(vertices_file_path::String)
    vertices = readdlm(vertices_file_path, ',', skipstart=1)
    normals = boundary_normals_calculator(vertices)
    return FreeFormBoundary(vertices, normals)
end


"""
    get_rectangle_boundary_vertices(farm_boundary::RectangleBoundary)

Returns the boundary vertices for a rectangular wind farm.

# Arguments
- `base::Float64`: length of the rectangle base
- `height::Float64`: length of the rectangle height
- `center::Array{Float64,1}`: position of the center of the rectangle (x,y)
"""
function get_rectangle_boundary_vertices(base, height, center)

    # calculate vertices
    vertices = [center[1]-base/2 center[2]-height/2;
                center[1]+base/2 center[2]-height/2;
                center[1]+base/2 center[2]+height/2;
                center[1]-base/2 center[2]+height/2]

    return vertices
end


"""
    get_rectangle_boundary_normals(farm_boundary::RectangleBoundary)

Returns the unit vectors normal to each boundary face for a rectangular wind farm.

# Arguments
- `base::Float64`: length of the rectangle base
- `height::Float64`: length of the rectangle height
- `center::Array{Float64,1}`: position of the center of the rectangle (x,y)
"""
function get_rectangle_boundary_normals(base, height, center)

    # calculate vertices
    vertices = get_rectangle_boundary_vertices(base, height, center)

    # calculate normals
    normals = boundary_normals_calculator(vertices)

    return normals
end


"""
    boundary_normals_calculator(boundary_vertices)

Outputs the unit vectors perpendicular to each boundary face, given the Cartesian coordinates for the shape's vertices.

# Arguments
- `boundary_vertices::Array{Float,1}` : n-by-2 array containing all the boundary vertices, counterclockwise
"""
function boundary_normals_calculator(boundary_vertices)
    # get number of vertices in shape
    nvertices = length(boundary_vertices[:, 1])
    # add the first vertex to the end of the array to form a closed loop
    boundary_vertices = [boundary_vertices; boundary_vertices[1,1] boundary_vertices[1,2]]
    # initialize array to hold boundary normals
    boundary_normals = zeros(nvertices, 2)
    # iterate over each boundary
    for i = 1:nvertices
        # create a vector normal to the boundary
        boundary_normals[i, :] = [-(boundary_vertices[i+1,2] - boundary_vertices[i,2]); boundary_vertices[i+1,1] - boundary_vertices[i,1]]
        # normalize the vector
        boundary_normals[i, :] = boundary_normals[i,:]/sqrt(sum(boundary_normals[i,:].^2))
    end

    return boundary_normals
end


"""
    generate_random_layout(layout_path, nturbines, rotor_diameter, boundary::CircleBoundary; min_spacing=1.0)

Generates a random circular wind farm layout and saves the coordinates to a text file.

# Arguments
- `layout_path::String`: path to the output text file that will contain the turbine layout coordinates
- `nturbines::Int64`: number of turbines in the farm
- `rotor_diameter::Float64`: rotor diameter of the turbines
- `boundary::CircleBoundary`: circular boundary structure containing the parameters for the boundary
- `min_spacing::Float64`: the minimum number of rotor diameters for the spacing between turbines
"""
function generate_random_layout(layout_path, nturbines, rotor_diameter, boundary::CircleBoundary; min_spacing=6.0)

    #################################################################################
    # SET UP
    #################################################################################

    # set bounding box
    xrange = [boundary.center[1] - boundary.radius; boundary.center[1] + boundary.radius]
    yrange = [boundary.center[2] - boundary.radius; boundary.center[2] + boundary.radius]

    # scale the minimum spacing
    min_spacing *= rotor_diameter


    #################################################################################
    # GENERATE RANDOM WIND FARM LAYOUT
    #################################################################################

    # create the boundary constraint function
    boundary_func(x, y) = ff.circle_boundary(boundary.center, boundary.radius, x, y)

    # generate the random layout
    turbine_x, turbine_y = generate_random_points_in_boundary(nturbines, boundary_func, xrange, yrange, min_spacing)


    #################################################################################
    # WRITE TURBINE COORDINATES TO A FILE
    #################################################################################

    # write turbine coordinates to a text file
    open(layout_path, "w") do io
        write(io, "# turbine_x, turbine_y\n")
        writedlm(io, [turbine_x turbine_y])
    end

end


"""
    generate_random_layout(layout_path, nturbines, rotor_diameter, boundary::Union{RectangleBoundary, FreeFormBoundary}; min_spacing=1.0)

Generates a random wind farm layout and saves the coordinates to a text file.

# Arguments
- `layout_path::String`: path to the output text file that will contain the turbine layout coordinates
- `nturbines::Int64`: number of turbines in the farm
- `rotor_diameter::Float64`: rotor diameter of the turbines
- `boundary::Union{RectangleBoundary, FreeFormBoundary}`: structure containing the parameters for the boundary
- `min_spacing::Float64`: the minimum number of rotor diameters for the spacing between turbines
"""
function generate_random_layout(layout_path, nturbines, rotor_diameter, boundary::Union{RectangleBoundary, FreeFormBoundary}; min_spacing=6.0)

    #################################################################################
    # SET UP
    #################################################################################

    # set bounding box
    xmin = minimum(boundary.vertices[:,1])
    xmax = maximum(boundary.vertices[:,1])
    ymin = minimum(boundary.vertices[:,2])
    ymax = maximum(boundary.vertices[:,2])
    xrange = [xmin; xmax]
    yrange = [ymin; ymax]

    # scale the minimum spacing
    min_spacing *= rotor_diameter


    #################################################################################
    # GENERATE RANDOM WIND FARM LAYOUT
    #################################################################################

    # create the boundary constraint function
    boundary_func(x, y) = ff.ray_trace_boundary(boundary.vertices, boundary.normals, x, y)

    # generate the random layout
    turbine_x, turbine_y = generate_random_points_in_boundary(nturbines, boundary_func, xrange, yrange, min_spacing)


    #################################################################################
    # WRITE TURBINE COORDINATES TO A FILE
    #################################################################################

    # write turbine coordinates to a text file
    open(layout_path, "w") do io
        write(io, "# turbine_x, turbine_y\n")
        writedlm(io, [turbine_x turbine_y])
    end

end


"""
    generate_random_points_in_boundary(n, boundary_func, xrange, yrange, min_spacing)

Generates n points randomly spaced inside of a boundary shape.

# Arguments
- `n::Int64`: number fo turbines in the farm
- `boundary_func::Function`: the boundary constraint function. Requires two arguments: (1) a vector of turbine x-coordinates and (2) a vector of turbine y-coordinates. Returns a vector of constraint values (one for each turbine).
- `xrange::Array{Float64,1}`: two-element vector with the min and max x-values for the random turbine x-coordinates
- `yrange::Array{Float64,1}`: two-element vector with the min and max values for the random turbine y-coordinates
- `min_spacing::Float64`: the minimum required distance between each turbine (not normalized)
"""
function generate_random_points_in_boundary(n, boundary_func, xrange, yrange, min_spacing)

    # initialize array to hold turbine coordinates
    locations = zeros(Float64, n, 2)
    # set counter to prevent while-loop from running forever
    count = 0
    for i = 1:n
        good_point = false
        while !good_point && count < 1e8
            # set good_point initially to true
            good_point = true   # this will be set to false if the generated point is not feasible
            # print out the number of loops every 1e6 loops
            if mod(count,1e6)==0 && count > 0
                println("Loop count: ", count, "\nLoosening the minimum spacing requirement.")
                min_spacing *= 0.75
            end
            # update the counter
            count += 1
            # generate random point in the bounding box
            locations[i,:] = 2*(rand(1,2) .- 0.5).*[xrange[2]-xrange[1] yrange[2]-yrange[1]]
            # calculate the distance(s) from the point to the boundary(s)
            distances_to_boundaries = boundary_func(locations[1:i,1], locations[1:i,2])
            # determine if the point is inside the wind farm boundary (+ distance is outside, - distance is inside)
            for j = 1:length(distances_to_boundaries)
                if distances_to_boundaries[j] > 0.0
                    good_point = false
                end
            end
            # determine if the point is far enough away from other points
            for turb = 1:i-1
                spacing = sqrt((locations[turb,1]-locations[i,1])^2 + (locations[turb,2]-locations[i,2])^2)
                if spacing < min_spacing
                    good_point = false
                end
            end
        end
    end
    x = locations[:,1]
    y = locations[:,2]

    return x, y
end

