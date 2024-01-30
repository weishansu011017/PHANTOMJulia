"""
Analysis Based on a specific ring with given r and z
"""

function spiral_analysis(data::phjlRawDataFrame,r::Float64, n_theta::Int, smoothed_kernal:: Function = M4_spline,h_mode::String="intep")
    """
    Make the data for the spiral analysis
    """

    @info "Start spiral analysis."
    kdtree2d = Generate_KDtree(data, 2)
    function wrap_surfdens(data::phjlRawDataFrame ,point::Array)
        return surface_density(data, point, smoothed_kernal,h_mode,"polar")
    end
    function wrap_quant(data::phjlRawDataFrame ,point::Array)
        return quantity_intepolate_2D(data, point,["tilt"],smoothed_kernal,h_mode,"polar")
    end
    if !(haskey(data.params, "h_mean"))
        add_mean_h(data)
    end

    Result_grid_dict = Dict{String,gridbackend}()
    Result_grid_dict["Sigma"] = generate_empty_grid([0.0],[2*pi],[n_theta])
    Result_grid_dict["tilt"] = generate_empty_grid([0.0],[2*pi],[n_theta])

    gridv = generate_coordinate_grid(Result_grid_dict["Sigma"])
    #Set up the truncated_radius
    truncated_radius = get_truncated_radius(data, 0.5, smoothed_kernal)
    for subarr in gridv
        pushfirst!(subarr, r)
    end

    for i in eachindex(gridv)
        target = gridv[i]
        kdtf_data = KDtree_filter(data, kdtree2d, target, truncated_radius,"polar")
        Result_grid_dict["Sigma"].grid[i] = wrap_surfdens(kdtf_data,target)
        Result_grid_dict["tilt"].grid[i] = wrap_quant(kdtf_data,target)["tilt"]
    end
    @info "End spiral analysis."
    return Result_grid_dict
end


export spiral_analysis