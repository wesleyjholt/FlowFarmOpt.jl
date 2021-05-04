abstract type AEPPlot end

struct BoxPlot{} <:AEPPlot
    # box plot for each group of data
end

struct BoxScatterPlot{} <: AEPPlot
    # box plot with overlayed scatter plot for each group of data
end

struct ConfidenceIntervalPlot <: AEPPlot
    # confidence interval plot
end

struct ConfidenceIntervalScatterPlot <: AEPPlot
    # confidence interval with scatter plot
end

struct SurfacePlot <: AEPPlot
    # surface plot (using Plots.jl package)
end

struct LinePlot <: AEPPlot
    # line plots for each layout
end


"""
    plot_initial_final_layout(x_initial, x_final, turbine_design, boundary::CircleBoundary; show_fig=true, save_fig=false, path_to_fig_file="layout.png")

Plots the initial and final layouts of circular wind farm.

# Arguments
- `x_initial::Array{Float64,1}`: initial turbine x and y coordinates, all concatenated into a single column vector
- `x_final::Array{Float64,1}`: final turbine x and y coordinates, all concatenated into a single column vector
- `turbine_design::TurbineType`: container holding all design parameters for a turbine type
- `boundary::CircleBoundary`: container holding parameters for a circular wind farm boundary
"""
function plot_initial_final_layout(x_initial, x_final, turbine_design, boundary::CircleBoundary; show_fig=true, save_fig=false, path_to_fig_file="layout.png")

    # add wind farm boundary to plot
    plt.gcf().gca().add_artist(plt.Circle((boundary.center[1],boundary.center[2]), boundary.radius, fill=false, color="black"))
    
    # set figure bounds
    axis("square")
    xlim(boundary.center[1] - boundary.radius*1.3, boundary.center[1] + boundary.radius*1.3)
    ylim(boundary.center[2] - boundary.radius*1.3, boundary.center[1] + boundary.radius*1.3)
    
    # plot turbines
    plot_turbines(x_initial, turbine_design.rotor_diameter, color="C1", alpha=0.1)
    plot_turbines(x_final, turbine_design.rotor_diameter)

    # save figure
    if save_fig
        savefig(path_to_fig_file, dpi=600)
    end

    # show figure
    if show_fig
        plt.show()
    end
end


function plot_initial_final_layout(x_initial, x_final, turbine_design, boundary::RectangleBoundary; show_fig=true, save_fig=false, path_to_fig_file="layout.png")

end


function plot_initial_final_layout(x_initial, x_final, turbine_design, boundary::FreeFormBoundary; show_fig=true, save_fig=false, path_to_fig_file="layout.png")

    # add wind farm boundary to plot
    boundary_vertices_closed = [boundary.vertices; reshape(boundary.vertices[1,:], 1, 2)]
    plt.gcf().gca().plot(boundary_vertices_closed[:,1], boundary_vertices_closed[:,2], color="black")
    
    # set figure bounds
    axis("square")
    xmin = minimum(boundary.vertices[:,1])
    xmax = maximum(boundary.vertices[:,1])
    ymin = minimum(boundary.vertices[:,2])
    ymax = maximum(boundary.vertices[:,2])
    xrange = xmax - xmin
    yrange = ymax - ymin
    plot_window_buffer = 0.1
    xlim(xmin - xrange*plot_window_buffer, xmax + xrange*plot_window_buffer)
    ylim(ymin - yrange*plot_window_buffer, ymax + yrange*plot_window_buffer)
    
    # plot turbines
    plot_turbines(x_initial, turbine_design.rotor_diameter, color="C1", alpha=0.1)
    plot_turbines(x_final, turbine_design.rotor_diameter)

    # save figure
    if save_fig
        savefig(path_to_fig_file, dpi=600)
    end

    # show figure
    if show_fig
        plt.show()
    end
end


function plot_overlayed_layouts()

end

"""
    plot_turbines(x, rotor_diameter; color="C0", fill=true, alpha=1.0, linestyle="-")

Plots turbines

# Arguments
- `x::Array{Float64,1}`: turbine x and y coordinates, all concatenated into a single column vector
- `rotor_diameter::Float64`: rotor diamter of the turbines
"""
function plot_turbines(x, rotor_diameter; color="C0", fill=true, alpha=1.0, linestyle="-")

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines] 
    turbine_y = x[nturbines+1:end]

    # add initial turbine location to plot
    for i = 1:length(turbine_x)
        plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=fill, color=color, alpha=alpha, linewidth=0.0))
    end

end


function plot_aeps(data_file_names::Array{String,2}, x_values, y_values, plot_type::SurfacePlot; fig_handle="", ax_handle="", show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="aep_boxplots.png", title="AEP Values for Optimized Layouts", xlabel="", ylabel="", upper_10percent=false)

    Plots.gr()
    # get the AEP values from the files
    aeps = get_aep_values_from_file_names(data_file_names)
    # filter out NaN values from AEP
    # aeps = filter_aeps(aeps_raw, layouts_must_match=true)
    # specify the number of groups
    n_x = length(x_values)
    n_y = length(y_values)
    # get x and y arrays
    X = x_values .* ones(n_x)
    Y = ones(n_y) .* y_values
    # get mean AEP value for each X-Y combination
    mean_aeps = zeros(n_ndirs_vec, n_nspeeds_vec)
    for i = 1:n_ndirs_vec
        for j = 1:n_nspeeds_vec
            if !upper_10percent
                mean_aeps[i,j] = mean(aeps[i,j])
            else
                mean_aeps[i,j] = calc_mean_upper10percent(aeps[i,j])
            end
        end
    end
    # normalize the mean AEP values
    norm_mean_aeps = mean_aeps ./ maximum(maximum.(aeps)) * 100.0
    println("Max AEP: $(maximum(maximum.(aeps)))")
    # reverse n_directions and n_speeds vectors to be in ascending order
    x_values_reversed = reverse(x_values)
    y_values_reversed = reverse(y_values)
    norm_mean_aeps_reversed = reverse(reverse(norm_mean_aeps, dims=1), dims=2)

    # interpolate mean AEP values onto a finer grid for plotting
    x_plot = range(minimum(x_values_reversed), stop=maximum(x_values_reversed), length=100)
    y_plot = range(minimum(y_values_reversed), stop=maximum(y_values_reversed), length=100)
    grid = RectangleGrid(x_values_reversed, y_values_reversed)
    z_plot = zeros(length(x_plot), length(y_plot))
    for i = 1:length(x_plot)
        for j = 1:length(y_plot)
            z_plot[j,i] = interpolate(grid, norm_mean_aeps_reversed, [y_plot[i], x_plot[j]])
        end
    end

    # # create heatmap
    # Plots.plot(x_plot, y_plot, z_plot, 
    #     st=:heatmap, 
    #     c=:blues, 
    #     background_color=:white, 
    #     foreground_color=:black, 
    #     top_margin = 10px,
    #     xlabel=xlabel,
    #     ylabel=ylabel,
    #     colorbar_title = "Normalized AEP (%)",
    #     framestyle=:box,
    #     clims=(minimum(norm_mean_aeps), maximum(norm_mean_aeps)),
    #     dpi=500
    #     )
    
    if !upper_10percent
        _clims = (98.6, 100.0)
        _colorbar_title = "Mean Normalized AEP (%)"
    else
        _clims = (99.0, 100.0) 
        _colorbar_title = "Mean of Upper 10% of Norm AEP (%)"
    end
    # create contour plot
    Plots.plot(x_plot, y_plot, (x,y) -> interpolate(grid, norm_mean_aeps_reversed, [x,y]), 
        st=:contour, 
        fill=true,
        c=:Blues_9, 
        background_color=:transparent, 
        foreground_color=:black, 
        foreground_color_border=:black,
        top_margin=25px,
        right_margin=10px,
        left_margin=10px,
        bottom_margin=10px,
        linewidth=0,
        xlabel=xlabel,
        ylabel=ylabel,
        xguidefontsize=15,
        yguidefontsize=15,
        colorbar_title = _colorbar_title,
        # colorbar_titlefonthalign = :right,
        colorbar_title_location = :top,
        framestyle=:box,
        clims=_clims,
        # cticks=98.6:.2:100.0,
        dpi=500
        )
    
    # save figure
    if save_fig
        mkpath(path_to_fig_directory)
        if !upper_10percent
            Plots.savefig(path_to_fig_directory * fig_file_name)
        else
            Plots.savefig(path_to_fig_directory * "upper-10percent-" * fig_file_name)
        end
    end
    # show figure
    if show_fig
        plt.show()
    end
end


"""
    plot_aeps(data_file_names, labels, plot_type::BoxPlot; show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="aep_boxplots.png", title="AEP Values for Optimized Layouts")
"""
function plot_aeps(data_file_names::Array{String,1}, x_values, plot_type::BoxPlot; fig_handle="", ax_handle="", show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="aep_boxplots.png", title="AEP Values for Optimized Layouts", xlabel="")

    # get figure and axes handles
    fig, ax = get_fig_ax_handles(fig_handle, ax_handle)
    # get the AEP values from the files
    aeps = get_aep_values_from_file_names(data_file_names)
    # create box plots
    create_box_plots(x_values, aeps, ax)
    # add plot labels
    add_aep_plot_labels(type=plot_type, title=title, xlabel=xlabel)
    # save figure
    if save_fig
        mkpath(path_to_fig_directory)
        savefig(path_to_fig_directory * fig_file_name , dpi=600)
    end
    # show figure
    if show_fig
        plt.show()
    end
end



"""
    plot_aeps(data_file_names, labels, plot_type::BoxScatterPlot; show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="aep_boxplots.png", title="AEP Values for Optimized Layouts")
"""
function plot_aeps(data_file_names::Array{String,1}, x_values, plot_type::BoxScatterPlot; fig_handle="", ax_handle="", show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="aep_boxplots_with_scatter.png", title="AEP Values for Optimized Layouts", xlabel="")

    # get figure and axes handles
    fig, ax = get_fig_ax_handles(fig_handle, ax_handle)
    # get the AEP values from the files
    aeps = get_aep_values_from_file_names(data_file_names)
    # create box plots
    create_box_plots(x_values, aeps, ax)
    # create scatter plots
    create_aep_scatter_plots(1:length(aeps), aeps, size=2)
    # add plot labels
    add_aep_plot_labels(type=plot_type, title=title, xlabel=xlabel)
    # save figure
    if save_fig
        mkpath(path_to_fig_directory)
        savefig(path_to_fig_directory * fig_file_name, dpi=600)
    end
    # show figure
    if show_fig
        plt.show()
    end
end



"""
    plot_aeps(data_file_names, labels, plot_type::ConfidenceIntervalPlot; show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="aep_boxplots.png", title="AEP Values for Optimized Layouts")
"""
function plot_aeps(data_file_names::Array{String,1}, x_values, plot_type::ConfidenceIntervalPlot; fig_handle="", ax_handle="", show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="aep_confidence_interval_plot.png", _title="AEP Values for Optimized Layouts", _xlabel="", upper_10percent=false)

    # get figure and axes handles
    fig, ax = get_fig_ax_handles(fig_handle, ax_handle)
    # get the AEP values from the files
    aeps = get_aep_values_from_file_names(data_file_names)
    # create mean line with confidence interval
    create_confidence_interval_plot(x_values, aeps, upper_10percent=upper_10percent)
    # add plot labels
    add_aep_plot_labels(type=plot_type, title=title, xlabel=xlabel)
    # save figure
    if save_fig
        mkpath(path_to_fig_directory)
        savefig(path_to_fig_directory * fig_file_name , dpi=600)
    end
    # show figure
    if show_fig
        plt.show()
    end
end


"""
    plot_aeps(data_file_names, labels, plot_type::ConfidenceIntervalScatterPlot; show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="aep_boxplots.png", title="AEP Values for Optimized Layouts")
"""
function plot_aeps(data_file_names::Array{String,1}, x_values, plot_type::ConfidenceIntervalScatterPlot; max_aep = 1.0, fig_handle="", ax_handle="", show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="aep_confidence_interval_plot_with_scatter.png", _title="AEP Values for Optimized Layouts", _xlabel="")

    # get figure and axes handles
    fig, ax = get_fig_ax_handles(fig_handle, ax_handle)
    # get the AEP values from the files
    aeps_raw = get_aep_values_from_file_names(data_file_names)
    # filter out NaN values
    aeps = filter_aeps(aeps_raw)/max_aep
    # create scatter plots
    create_aep_scatter_plots(x_values, aeps)
    # create mean line with confidence interval
    create_confidence_interval_plot(x_values, aeps)
    # add plot labels
    title(_title)
    xlabel(_xlabel)
    if max_aep == 1.0
        ylabel("AEP (MWh)")
    else
        ylabel("Normalized AEP")
        ylim([0.979, 1.0])
        yticks([0.98, 0.985, 0.99, 0.995, 1.0])
    end
    # save figure
    if save_fig
        mkpath(path_to_fig_directory)
        savefig(path_to_fig_directory * fig_file_name , dpi=600)
    end
    # show figure
    if show_fig
        plt.show()
    end
end


function plot_aeps(data_file_names::Array{String,1}, x_values, plot_type::LinePlot; fig_handle="", ax_handle="", show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="aep_line_plot.png", title="AEP Values for Optimized Layouts", xlabel="")

    # get figure and axes handles
    fig, ax = get_fig_ax_handles(fig_handle, ax_handle)
    # get the AEP values from the files
    aeps = get_aep_values_from_file_names(data_file_names)
    # create line plots
    create_line_plot(x_values, aeps)
    # add plot labels
    add_aep_plot_labels(type=plot_type, title=title, xlabel=xlabel)
    # save figure
    if save_fig
        mkpath(path_to_fig_directory)
        savefig(path_to_fig_directory * fig_file_name , dpi=600)
    end
    # show figure
    if show_fig
        plt.show()
    end
end


function create_aep_surface(x_values, y_values, aeps)

    # specify the number of groups
    n_x = length(x_values)
    n_y = length(y_values)
    # get x and y arrays
    X = x_values .* ones(n_x)
    Y = ones(n_y) .* y_values
    # get mean AEP value for each X-Y combination
    mean_aeps = zeros(n_ndirs_vec, n_nspeeds_vec)
    for i = 1:n_ndirs_vec
        for j = 1:n_nspeeds_vec
            mean_aeps[i,j] = mean(aeps[i,j])
        end
    end
    # create surface plot
    Plots.plot(x_values, y_values, mean_aeps, st=:surface, c=:blues, camera=(90,90), zaxis=false)
    # ax = plt.gca()
    # ax.view_init(90, 90)
    # ax.w_zaxis.line.set_lw(0.)
    # ax.set_zticks([])
    # plt.colorbar(colormap)
    # plt.gcf().colorbar(surf, shrink=0.5, aspect=5)
end


function create_box_plots(x_values, aeps, ax)
    flierprops = Dict{String, Any}()
    flierprops["markersize"] = 3
    flierprops["markeredgewidth"] = 0.5
    plt.boxplot(aeps, whis=(5, 95), flierprops=flierprops, zorder=1)
    if isa(x_values, Array{Float64,1})
        ax.set_xticklabels(round.(x_values, digits=1))
    else
        ax.set_xticklabels(x_values)
    end
    ax.xaxis.grid(true, alpha=0.2, linestyle="--")
end



function create_aep_scatter_plots(x_values, aeps; alpha=0.3, size=2, jitter_standard_dev=0.05)
    for i = 1:length(aeps)
        x = zeros(length(aeps[i])) .+ x_values[i] .+ randn(length(aeps[i]))*jitter_standard_dev
        plt.scatter(x, aeps[i], color="C0", alpha=alpha, s=size, zorder=10)
    end
end

function create_confidence_interval_plot(x_values, aeps; confidence_interval_width=1.0, upper_10percent=false)
    if !upper_10percent
        aeps_means = mean.(aeps)
        aeps_std = std.(aeps)
    else
        aeps_means = calc_mean_upper10percent.(aeps)
        aeps_std = calc_std_upper10percent.(aeps)
    end
    aeps_lower_CI = aeps_means - aeps_std*confidence_interval_width
    aeps_upper_CI = aeps_means + aeps_std*confidence_interval_width
    plt.plot(x_values, aeps_means, color="black", alpha=0.5, linewidth=1)
    aeps_lower_CI = reshape(aeps_lower_CI, length(aeps_lower_CI))
    aeps_upper_CI = reshape(aeps_upper_CI, length(aeps_upper_CI))
    plt.fill_between(x_values, aeps_lower_CI, aeps_upper_CI, color="black", alpha=0.1, linewidth=0.0)
end

function create_line_plot(x_values, aeps)
    for i = 1:length(aeps[1])
        plt.plot(x_values, aeps[i], color="C0", alpha=0.5)
    end
end

function add_aep_plot_labels(; type=type::Union{BoxPlot, BoxScatterPlot, ConfidenceIntervalPlot, ConfidenceIntervalScatterPlot}, title="", x_label="", y_label="AEP (MWh)")
    title(title)
    xlabel(x_label)
    ylabel(y_label)
end

function add_aep_plot_labels(; type=type::SurfacePlot, title="", xlabel="", ylabel="")
    Plots.plot(title=title)
    # Plots.xlabel(xlabel)
    # Plots.ylabel(ylabel)
end


function get_fig_ax_handles(fig_handle, ax_handle)

    if fig_handle==""
        fig = plt.figure()
    else
        fig = fig_handle
    end
    
    if ax_handle==""
        ax = fig.add_axes([0.15, 0.1, 0.8, 0.8])
    else
        ax = ax_handle
    end

    return fig, ax
end



function plot_montecarlo(data_file_names::Array{String,2}, max_aep, sample_sizes, ndirs_vec, nspeeds_vec; fig_handle="", ax_handle="", show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="montecarlo_plot.png", title="Monte Carlo Simulation", xlabel="")

    # get the AEP values from the files
    aeps = get_aep_values_from_file_names(data_file_names)

    # throw out AEPs for layouts that have at least one NaN
    # aeps_filtered = filter_aeps(aeps_raw)

    # normalize AEP values (assuming aeps_filtered[end,end] contains results for 360 directions and 25 wind speeds)
    # aeps_normalized = normalize_aeps(aeps_filtered)
    # println(aeps_normalized)

    # get number of directions and speeds
    n_ndirs = length(ndirs_vec)
    n_nspeeds = length(nspeeds_vec)

    # initialize convergence plots
    figure()
    subplot(211)
    plt.title("Monte Carlo Convergence Plots")
    plt.ylabel("Mean AEP\n(MWh)")

    subplot(212)
    plt.xlabel("Sample Size")
    plt.ylabel("Standard Deviation AEP\n(MWh)")

    # iterate through each speed and direction
    legend_labels = fill("", Int(n_ndirs*n_nspeeds))
    for i = 1:n_ndirs
        for j = 1:n_nspeeds
            # group into samples and calculate means and standard deviations
            n_sample_sizes = length(sample_sizes[i,j])
            mean_aeps = zeros(n_sample_sizes)
            std_aeps = zeros(n_sample_sizes)
            index_1 = 1
            index_2 = 0
            for k = 1:n_sample_sizes
                index_2 += Int(sample_sizes[i,j][k])
                mean_aeps[k] = mean(aeps[i,j][index_1:index_2])
                std_aeps[k] = std(aeps[i,j][index_1:index_2])
                index_1 = deepcopy(index_2)
            end
            mean_aeps_normalized = mean_aeps/max_aep
            subplot(211)
            scatter(sample_sizes[i,j], mean_aeps_normalized, color="C0")
            subplot(212)
            scatter(sample_sizes[i,j], std_aeps, color="C0")
            legend_labels[(i-1)*n_nspeeds + j] = "$(ndirs_vec[i]) dirs, $(nspeeds_vec[j]) speeds"
        end
    end

    # add legend
    plt.legend(legend_labels)

    # save figure
    if save_fig
        savefig(path_to_fig_directory * fig_file_name, dpi=600)
    end

    # show figure
    if show_fig
        plt.show()
    end
end



function plot_montecarlo_old(data_file_name::String, max_aep, sample_sizes, ndirs_vec, nspeeds_vec; fig_handle="", ax_handle="", show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="montecarlo_plot.png", title="", xlabel="")

    # get the AEP values from the files
    aeps = get_aep_values_from_file_names(data_file_name)

    # initialize convergence plots
    fig = figure()
    ax1 = fig.add_axes([0.18, 0.55, 0.8, 0.35])
    plt.title(title, fontsize=16)
    plt.ylabel("Mean Normalized AEP\n(MWh)", fontsize=14)

    ax2 = fig.add_axes([0.18, 0.1, 0.8, 0.35])
    plt.xlabel("Sample Size", fontsize=14)
    plt.ylabel("95th Percentile AEP\n(MWh)", fontsize=14)

    # group into samples and calculate means and standard deviations
    n_sample_sizes = length(sample_sizes)
    mean_aeps = zeros(n_sample_sizes)
    mean_upper10percent_aeps = zeros(n_sample_sizes)
    index_1 = 1
    index_2 = 0
    for k = 1:n_sample_sizes
        index_2 += Int(sample_sizes[k])
        mean_aeps[k] = mean(aeps[index_1:index_2])
        mean_upper10percent_aeps[k] = calc_mean_upper10percent(aeps[index_1:index_2])
        index_1 = deepcopy(index_2)
    end
    mean_aeps_normalized = mean_aeps/max_aep
    mean_upper10percent_aeps_normalized = mean_upper10percent_aeps/max_aep
    
    ax1.scatter(sample_sizes, mean_aeps_normalized, color="C0")
    # yticks([floor(minimum(mean_aeps_normalized), digits=4):0.0001:ceil(maximum(mean_aeps_normalized), digits=4)])
    
    ax2.scatter(sample_sizes, mean_upper10percent_aeps_normalized, color="C0")
    yticks([floor(minimum(mean_upper10percent_aeps), digits=4):0.0001:ceil(maximum(mean_upper10percent_aeps), digits=4)])
    # ylim([0.9979, 0.9984])

    # save figure
    if save_fig
        savefig(path_to_fig_directory * fig_file_name, dpi=600)
    end

    # show figure
    if show_fig
        plt.show()
    end
end


function plot_montecarlo(data_file_name::String, max_aep, sample_sizes, ndirs_vec, nspeeds_vec; fig_handle="", ax_handle="", show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="montecarlo_plot.png", title="", xlabel="")

    # get the AEP values from the files
    aeps = get_aep_values_from_file_names(data_file_name)

    # initialize scatter plot
    fig = figure()
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8])
    plt.title(title, fontsize=16)
    plt.xlabel("Sample Size", fontsize=14)
    plt.ylabel("Mean Normalized AEP", fontsize=14)

    # group into samples and calculate means and standard deviations
    n_sample_sizes = length(sample_sizes)
    mean_aeps = zeros(n_sample_sizes)
    mean_upper10percent_aeps = zeros(n_sample_sizes)
    index_1 = 1
    index_2 = 0
    for k = 1:n_sample_sizes
        index_2 += Int(sample_sizes[k])
        mean_aeps[k] = mean(aeps[index_1:index_2])
        mean_upper10percent_aeps[k] = calc_mean_upper10percent(aeps[index_1:index_2])
        index_1 = deepcopy(index_2)
    end
    mean_aeps_normalized = mean_aeps/max_aep
    mean_upper10percent_aeps_normalized = mean_upper10percent_aeps/max_aep
    
    ax.scatter(sample_sizes, mean_aeps_normalized, color="C0", alpha=0.5, edgecolors="none")
    # yticks([floor(minimum(mean_aeps_normalized), digits=4):0.0001:ceil(maximum(mean_aeps_normalized), digits=4)])
    
    ax.scatter(sample_sizes, mean_upper10percent_aeps_normalized, color="C1", alpha=0.5, edgecolors="none")
    # yticks([floor(minimum(mean_upper10percent_aeps), digits=4):0.0001:ceil(maximum(mean_upper10percent_aeps), digits=4)])
    ylim([0.997, 1.0])

    legend(["Mean of entire sample", "Mean of upper 10% of sample"])

    # save figure
    if save_fig
        savefig(path_to_fig_directory * fig_file_name, dpi=600)
    end

    # show figure
    if show_fig
        plt.show()
    end
end


function plot_initial_final_aep_scatter(initial_data_file_name::String, final_data_file_name::String, max_aep; fig_handle="", ax_handle="", show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="initial_final_aep_plot.png", title="", xlabel="")

    # get the AEP values from the files
    initial_aeps = get_aep_values_from_file_names(initial_data_file_name)
    final_aeps = get_aep_values_from_file_names(final_data_file_name)

    # initialize scatter plot
    fig = figure()
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8])
    plt.title(title, fontsize=16)
    plt.xlabel("Initial Layout Normalized AEP", fontsize=14)
    plt.ylabel("Final Layout Normalized AEP", fontsize=14)
    plt.xlim([0.93,1.0])
    plt.ylim([0.994,1.0])
    println(length(initial_aeps))
    println(length(final_aeps))
    # plot the data
    ax.scatter(initial_aeps ./ max_aep, final_aeps ./ max_aep, color="C0", s = 8, alpha=.3, edgecolors="none")

    # save figure
    if save_fig
        savefig(path_to_fig_directory * fig_file_name, dpi=600)
    end

    # show figure
    if show_fig
        plt.show()
    end

end


function plot_initial_final_aep_scatter_categorical(initial_data_file_name::String, final_data_file_name::String; fig_handle="", ax_handle="", show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="initial_final_aep_plot.png", title="", xlabel="")

    # get the AEP values from the files
    initial_aeps = get_aep_values_from_file_names(initial_data_file_name)
    final_aeps = get_aep_values_from_file_names(final_data_file_name)

    # initialize scatter plot
    fig = figure()
    ax = fig.add_axes([0.15, 0.1, 0.8, 0.8])
    plt.title(title, fontsize=16)
    plt.xlabel("Initial Layout Percentile", fontsize=14)
    plt.ylabel("Final Layout Percentile", fontsize=14)

    # get percentile values
    x_percentile = (sortperm(initial_aeps) .- 1)/length(initial_aeps) * 100
    y_percentile = (sortperm(final_aeps) .- 1)/length(initial_aeps) * 100

    # plot the data
    ax.scatter(x_percentile, y_percentile, color="C0", s = 8, alpha=.3, edgecolors="none")

    # save figure
    if save_fig
        savefig(path_to_fig_directory * fig_file_name, dpi=600)
    end

    # show figure
    if show_fig
        plt.show()
    end

end


# function plot_montecarlo(data_file_names::Array{String,2}, sample_sizes, ndirs_vec, nspeeds_vec; fig_handle="", ax_handle="", show_fig=false, save_fig=true, path_to_fig_directory="", fig_file_name="montecarlo_plot.png", title="Monte Carlo Simulation", xlabel="")

#     # get the AEP values from the files
#     aeps = get_aep_values_from_file_names(data_file_names)

#     # get number of directions and speeds
#     n_ndirs = length(ndirs_vec)
#     n_nspeeds = length(nspeeds_vec)

#     # initialize convergence plots
#     figure()
#     subplot(211)
#     plt.title("Monte Carlo Convergence Plots")
#     plt.ylabel("Mean AEP\n(MWh)")

#     subplot(212)
#     plt.xlabel("Sample Size")
#     plt.ylabel("Standard Deviation AEP\n(MWh)")

#     # iterate through each speed and direction
#     legend_labels = fill("", Int(n_ndirs*n_nspeeds))
#     for i = [1]#:n_ndirs
#         for j = [1]#:n_nspeeds
#             # group into samples and calculate means and standard deviations
#             n_sample_sizes = length(sample_sizes[i,j])
#             mean_aeps = zeros(n_sample_sizes)
#             std_aeps = zeros(n_sample_sizes)
#             index_1 = 1
#             index_2 = 0
#             for k = 1:n_sample_sizes
#                 index_2 += Int(sample_sizes[i,j][k])
#                 mean_aeps[k] = mean(aeps[i,j][index_1:index_2])
#                 std_aeps[k] = std(aeps[i,j][index_1:index_2])
#                 index_1 = deepcopy(index_2)
#             end
#             subplot(211)
#             scatter(sample_sizes[i,j], mean_aeps, color="C0")
#             subplot(212)
#             scatter(sample_sizes[i,j], std_aeps, color="C0")
#             legend_labels[(i-1)*n_nspeeds + j] = "$(ndirs_vec[i]) dirs, $(nspeeds_vec[j]) speeds"
#         end
#     end

#     # add legend
#     plt.legend(legend_labels)

#     # save figure
#     if save_fig
#         savefig(path_to_fig_directory * fig_file_name, dpi=600)
#     end

#     # show figure
#     if show_fig
#         plt.show()
#     end
# end

function pairs(data)
    (nobs, nvars) = size(data)
    (fig, ax) = subplots(nvars, nvars, figsize=(8,8))
    subplots_adjust(hspace=0.05, wspace=0.05)

    # Plot data
    for i = 1:nvars
    for j = 1:nvars
        if i != j
            ax[i,j][:plot](data[:,j],data[:,i],"ob",mfc="none")
        else
            ax[i,j][:hist](data[:,i])
        end
        ax[i,j][:xaxis][:set_visible](false)
        ax[i,j][:yaxis][:set_visible](false)
    end
    end

    # Set tick positions
    for i = 1:nvars
    ax[i,1][:yaxis][:set_ticks_position]("left")
    ax[i,end][:yaxis][:set_ticks_position]("right")
    ax[1,i][:xaxis][:set_ticks_position]("top")
    ax[end,i][:xaxis][:set_ticks_position]("bottom")
    end

    # Turn ticks on
    cc = repmat([nvars, 1],int(ceil(nvars/2)),1)
    for i = 1:nvars
    ax[i,cc[i]][:yaxis][:set_visible](true)
    ax[cc[i],i][:xaxis][:set_visible](true)
    end
end

