function density_analysis(data::phjlRawDataFrame,r_params::Tuple{Float64,Float64,Int}, z::Float64, n_theta :: Int, smoothed_kernal:: Function = M4_spline,h_mode::String="intep")
    """
    Calculate the density for a whole disc with a specific z
    """

    @info "Start density analysis."
    function wrap_dens(point)
        return density(data, point, smoothed_kernal,h_mode,"cylin")
    end
    imin = [r_params[1], 0.0]
    imax = [r_params[2], 2*pi]
    in = [r_params[3],n_theta]
    rho_grid = disc_2d_grid_generator(imin,imax,in)
    storage = rho_grid.grid
    gridv = generate_coordinate_grid(rho_grid)
    push!.(gridv,z)
    @threads for i in eachindex(gridv)
        storage[i] = wrap_dens(gridv[i])
    end
    @info "End density analysis."
    return rho_grid
end

function surface_density_analysis(data::phjlRawDataFrame,r_params::Tuple{Float64,Float64,Int}, n_theta :: Int, smoothed_kernal:: Function = M4_spline,h_mode::String="intep")
    """
    Calculate the density for a whole disc with z = 0
    """

    @info "Start surface density analysis."
    kdtree2d = Generate_KDtree(data, 2)
    function wrap_surfdens(data::phjlRawDataFrame ,point::Array)
        return surface_density(data, point, smoothed_kernal,h_mode,"polar")
    end
    if !(haskey(data.params, "h_mean"))
        add_mean_h(data)
    end
    imin = [r_params[1], 0.0]
    imax = [r_params[2], 2*pi]
    in = [r_params[3],n_theta]
    Sigma_grid = disc_2d_grid_generator(imin,imax,in)
    storage = Sigma_grid.grid
    gridv = generate_coordinate_grid(Sigma_grid)
    #Set up the truncated_radius
    truncated_radius = get_truncated_radius(data, 0.5, smoothed_kernal)
    @threads for i in eachindex(gridv)
        target = gridv[i]
        kdtf_data2d = KDtree_filter(data, kdtree2d, target, truncated_radius,"polar")
        storage[i] = wrap_surfdens(kdtf_data2d,target)
    end
    @info "End surface density analysis."
    return Sigma_grid
end

function midH_analysis(data::phjlRawDataFrame,r_params::Tuple{Float64,Float64,Int} ,n_theta :: Int, mid_frac::Float64 = 1.0, smoothed_kernal:: Function = M4_spline,h_mode::String="intep")
    """
    Calculate the midplane scale height.
    """

    @info "Start midplane scale height analysis."
    kdtree2d = Generate_KDtree(data, 2)
    kdtree3d = Generate_KDtree(data, 3)
    function wrap_dens(data::phjlRawDataFrame ,point::Array)
        return density(data, point, smoothed_kernal,h_mode,"cylin")
    end
    function wrap_surfdens(data::phjlRawDataFrame ,point::Array)
        return surface_density(data, point, smoothed_kernal,h_mode,"polar")
    end
    if !(haskey(data.params, "h_mean"))
        add_mean_h(data)
    end
    r2pi = 1/sqrt(2*pi)
    imin = [r_params[1], 0.0]
    imax = [r_params[2], 2*pi]
    in = [r_params[3],n_theta]

    midH_grid = disc_2d_grid_generator(imin,imax,in)
    Sigma_grid = deepcopy(midH_grid)
    gridv = generate_coordinate_grid(midH_grid)

    midH_storage = midH_grid.grid
    Sigma_storage = Sigma_grid.grid
    midH_array = zeros(Float64,in[1])
    #Set up the truncated_radius
    truncated_radius = get_truncated_radius(data, 0.5, smoothed_kernal)
    @threads for i in eachindex(gridv)
        target = gridv[i]
        target3d = [target..., 0.0]
        # Searching the Neighbors with a truncate_radius
        kdtf_data2d = KDtree_filter(data, kdtree2d, target, truncated_radius,"polar")
        kdtf_data3d = KDtree_filter(data, kdtree3d, target3d, truncated_radius,"polar")

        Sigmai = wrap_surfdens(kdtf_data2d,target)
        rho0i = wrap_dens(kdtf_data3d,target3d)
        if (Sigmai == 0.0) || (rho0i == 0.0)
            continue
        else
            Hi = Sigmai/(r2pi*rho0i)
            Sigma_storage[i] = Sigmai
            midH_storage[i] = Hi * mid_frac
        end
    end
    for i in eachindex(midH_array)
        midH_array[i] = mean(midH_storage[i,:])
    end
    @info "End midplane scale height analysis."
    return midH_array, Sigma_grid
end

function pitch_analysis(data::phjlRawDataFrame,r_params::Tuple{Float64,Float64,Int} ,n_theta :: Int, midH_array::Vector, z_saperate::Int =5, smoothed_kernal:: Function = M4_spline,h_mode::String="intep")
    """
    Make the data for pitch analysis.
    """
    @info "Start pitch analysis."
    kdtree3d = Generate_KDtree(data, 3)
    column_names = ["vr","vphi"]
    function wrap_dens(data::phjlRawDataFrame,point::Array)
        return density(data, point, smoothed_kernal,h_mode,"cylin")
    end
    function wrap_quant(data::phjlRawDataFrame,point::Array)
        return quantity_intepolate(data, point, column_names, smoothed_kernal,h_mode,"cylin")
    end

    if !(hasproperty(data.dfdata, "vr"))
        if (data.params["Origin_located"] == "COM")
            error("IntepolateError: Wrong origin located!")
        end
        add_cylindrical(data)
    end
    if !(haskey(data.params, "h_mean"))
        add_mean_h(data)
    end

    #Preparing the parameters for grid generator
    r2pi = 1/sqrt(2*pi)
    imin = [r_params[1], 0.0]
    imax = [r_params[2], 2*pi]
    in = [r_params[3],n_theta]

    # Generate a grid for midH
    midH_grid = disc_2d_grid_generator(imin,imax,in)
    midH_storage = midH_grid.grid
    for i in eachindex(midH_array)
        midH_storage[i,:] .= midH_array[i]
    end

    # Generate the gird for each analysis
    Result_grid_dict = Dict{String, gridbackend}()
    Result_grid_dict["rho_m"] = disc_2d_grid_generator(imin,imax,in)
    for column_name in column_names
        Result_grid_dict[column_name] = deepcopy(Result_grid_dict["rho_m"])
    end

    # Generate a coordinate grid
    gridv = generate_coordinate_grid(midH_grid)

    # Generate a storage array to store data
    Result_storage_dict = Dict{String, Array}()
    Result_storage_dict["rho_m"] = Result_grid_dict["rho_m"].grid
    for column_name in column_names
        Result_storage_dict[column_name] = Result_grid_dict[column_name].grid
    end

    truncated_radius = get_truncated_radius(data, 0.5, smoothed_kernal)
    @threads for i in eachindex(gridv)
        target = gridv[i]
        midHi = midH_storage[i]
        if (midHi == 0.0)
            target3d = [target..., 0.0]
            # Searching the Neighbors with a truncate_radius
            kdtf_data = KDtree_filter(data, kdtree3d, target3d, truncated_radius,"polar")
            # rhom
            Result_storage_dict["rho_m"][i] = wrap_dens(kdtf_data,target3d)
            # other quantity
            buffer_dict = wrap_quant(kdtf_data,target3d)
            for column_name in column_names
                Result_storage_dict[column_name][i] = buffer_dict[column_name]
            end
        else
            intepo_result = Dict{String, Vector}()
            for column_name in column_names
                intepo_result[column_name] = zeros(Float64,z_saperate)
            end
            buffer_rho = zeros(Float64,z_saperate)
            midthick = 2*midHi
            mid_array = LinRange(-midHi,midHi,z_saperate)
            for j in eachindex(mid_array)
                target3d = [target..., mid_array[j]]
                # Searching the Neighbors with a truncate_radius
                kdtf_data = KDtree_filter(data, kdtree3d, target3d, truncated_radius,"polar")
                # rhom
                buffer_rho[j] = wrap_dens(kdtf_data,target3d)
                # other quantity
                buffer_dict = wrap_quant(kdtf_data,target3d)
                for column_name in column_names
                    intepo_result[column_name][j] = buffer_dict[column_name]
                end
            end
            Result_storage_dict["rho_m"][i] = _Integral_1d(mid_array,buffer_rho,[-midHi,midHi])/midthick
            for column_name in column_names
                Result_storage_dict[column_name][i] = mean(intepo_result[column_name])
            end
        end
    end
    @info "End pitch analysis."
    return Result_grid_dict
end

export density_analysis, surface_density_analysis, midH_analysis, mid_density_analysis, pitch_analysis