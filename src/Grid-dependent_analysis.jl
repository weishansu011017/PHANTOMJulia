function density_analysis(data::phjlRawDataFrame,r_params::Tuple{Float64,Float64,Int}, z::Float64, n_theta :: Int, smoothed_kernal:: Function = M4_spline,h_mode::String="intep")
    """
    Calculate the density for a whole disc with a specific z in the cylindrical coordinate
    --------
    Step 1. Generate a 3d KDTree to speed up the searching of particles

    Step 2. Generate a grid for analysis. The grid is base on a struct gridbackend which is defined in the "grid.jl" file

    Step 3. Calculate the result for each point by using the SPH intepolation. The calculation for each point is defined in "physical_quantity.jl"
    --------
    """
    @info "Start density analysis."
    kdtree3d = Generate_KDtree(data, 3)
    function wrap_dens(data::phjlRawDataFrame ,point::Array)
        return density(data, point, smoothed_kernal,h_mode,"polar")
    end
    imin = [r_params[1], 0.0]
    imax = [r_params[2], 2*pi]
    in = [r_params[3],n_theta]
    rho_grid = disc_2d_grid_generator(imin,imax,in)
    storage = rho_grid.grid
    gridv = generate_coordinate_grid(rho_grid)
    push!.(gridv,z)
    @threads for i in eachindex(gridv)
        target = gridv[i]
        kdtf_data3d = KDtree_filter(data, kdtree3d, target, truncated_radius,"polar")
        storage[i] = wrap_dens(kdtf_data3d,gridv[i])
    end
    @info "End density analysis."
    return rho_grid
end

function surface_density_analysis(data::phjlRawDataFrame,r_params::Tuple{Float64,Float64,Int}, n_theta :: Int, smoothed_kernal:: Function = M4_spline,h_mode::String="intep")
    """
    Calculate the surface density for a whole disc 
    --------
    Step 1. Generate a 2d KDTree to speed up the searching of particles

    Step 2. Generate a grid for analysis. The grid is base on a struct gridbackend which is defined in the "grid.jl" file

    Step 3. Calculate the result for each point by using the SPH intepolation. The calculation for each point is defined in "physical_quantity.jl"
    --------
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
    --------
    Step 1. Generate a 2d and 3d KDTree to speed up the searching of particles

    Step 2. Generate a grid for analysis. The grid is base on a struct gridbackend which is defined in the "grid.jl" file

    Step 3. Calculate the density and surface density for each point by using the SPH intepolation. The calculation for each point is defined in "physical_quantity.jl"

    Step 4. Calculate the midplane scale height by the formula H_mid = mid_frac*(Sigma_g/ sqrt(2π)* ρ_g(z=0))
    --------
    """
    @info "Start midplane scale height analysis."
    kdtree2d = Generate_KDtree(data, 2)
    kdtree3d = Generate_KDtree(data, 3)
    function wrap_dens(data::phjlRawDataFrame ,point::Array)
        return density(data, point, smoothed_kernal,h_mode,"polar")
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
    gridv = generate_coordinate_grid(midH_grid)

    midH_storage = midH_grid.grid
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
            midH_storage[i] = Hi * mid_frac
        end
    end
    for i in eachindex(midH_array)
        midH_array[i] = mean(midH_storage[i,:])
    end
    @info "End midplane scale height analysis."
    return midH_array
end

function pitch_analysis(data::phjlRawDataFrame,r_params::Tuple{Float64,Float64,Int} ,n_theta :: Int, midH_array::Vector, z_saperate::Int,column_names::Vector{String}, smoothed_kernal:: Function = M4_spline,h_mode::String="intep")
    """
    Make the data for pitch analysis.
    --------
    Step 1. Generate a 3d KDTree to speed up the searching of particles

    Step 2. Generate a grid for analysis. The grid is base on a struct gridbackend which is defined in the "grid.jl" file

    Step 3. Calculate the result for each point by using the SPH intepolation. The calculation for each point is defined in "physical_quantity.jl"
    --------
    """
    @info "Start pitch analysis."
    kdtree3d = Generate_KDtree(data, 3)
    function wrap_dens(data::phjlRawDataFrame,point::Array)
        return density(data, point, smoothed_kernal,h_mode,"polar")
    end
    function wrap_quant(data::phjlRawDataFrame,point::Array)
        return quantity_intepolate(data, point, column_names, smoothed_kernal,h_mode,"polar")
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
    Result_grid_dict["rho"] = disc_2d_grid_generator(imin,imax,in)
    for column_name in column_names
        Result_grid_dict[column_name] = deepcopy(Result_grid_dict["rho"])
    end

    # Generate a coordinate grid
    gridv = generate_coordinate_grid(midH_grid)

    # Generate a storage array to store data
    Result_storage_dict = Dict{String, Array}()
    Result_storage_dict["rho"] = Result_grid_dict["rho"].grid
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
            Result_storage_dict["rho"][i] = wrap_dens(kdtf_data,target3d)
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
            Result_storage_dict["rho"][i] = _Integral_1d(mid_array,buffer_rho,[-midHi,midHi])/midthick
            for column_name in column_names
                Result_storage_dict[column_name][i] = _Integral_1d(mid_array,intepo_result[column_name],[-midHi,midHi])/midthick #_Integral_1d(mid_array,intepo_result[column_name],[-midHi,midHi])/midthick or mean(intepo_result[column_name]) 
            end
        end
    end
    @info "End pitch analysis."
    return Result_grid_dict
end

function pitch_other_analysis(data::phjlRawDataFrame,r_params::Tuple{Float64,Float64,Int}, n_theta :: Int, column_names::Vector{String}, smoothed_kernal:: Function = M4_spline,h_mode::String="intep")
    """
    Make the data for pitch analysis but no need to considered the midplane
    --------
    Step 1. Generate a 2d KDTree to speed up the searching of particles

    Step 2. Generate a grid for analysis. The grid is base on a struct gridbackend which is defined in the "grid.jl" file

    Step 3. Calculate the result for each point by using the SPH intepolation. The calculation for each point is defined in "physical_quantity.jl"
    --------
    """
    @info "Start 2d intepolate analysis."
    kdtree2d = Generate_KDtree(data, 2)
    function wrap_surfdens(data::phjlRawDataFrame ,point::Array)
        return surface_density(data, point, smoothed_kernal,h_mode,"polar")
    end
    function wrap_grad_surfdens(data::phjlRawDataFrame, point::Array)
        return gradient_surface_density(data, point, smoothed_kernal,h_mode,"polar")
    end
    function wrap_quant2d(data::phjlRawDataFrame,point::Array, Sigmai::Float64)
        return quantity_intepolate_2D(data, point,Sigmai, column_names, smoothed_kernal,h_mode,"polar")
    end
    if !(haskey(data.params, "h_mean"))
        add_mean_h(data)
    end
    if !(hasproperty(data.dfdata, "e"))
        error("ColumnMissingError: Please add the eccentricity column first.")
    end
    imin = [r_params[1], 0.0]
    imax = [r_params[2], 2*pi]
    in = [r_params[3],n_theta]
    # Generate the gird for each analysis
    Result_grid_dict = Dict{String, gridbackend}()
    Result_grid_dict["Sigma"] = disc_2d_grid_generator(imin,imax,in)
    for column_name in column_names
        Result_grid_dict[column_name] = deepcopy(Result_grid_dict["Sigma"])
    end
    Result_grid_dict["∇Sigmax"] = deepcopy(Result_grid_dict["Sigma"])
    Result_grid_dict["∇Sigmay"] = deepcopy(Result_grid_dict["Sigma"])

    # Generate a coordinate grid
    gridv = generate_coordinate_grid(Result_grid_dict["Sigma"])

    # Generate a storage array to store data
    Result_storage_dict = Dict{String, Array}()
    Result_storage_dict["Sigma"] = Result_grid_dict["Sigma"].grid
    for column_name in column_names
        Result_storage_dict[column_name] = Result_grid_dict[column_name].grid
    end
    Result_storage_dict["∇Sigmax"] = Result_grid_dict["∇Sigmax"].grid
    Result_storage_dict["∇Sigmay"] = Result_grid_dict["∇Sigmay"].grid

    #Set up the truncated_radius
    truncated_radius = get_truncated_radius(data, 0.5, smoothed_kernal)
    @threads for i in eachindex(gridv)
        target = gridv[i]
        kdtf_data2d = KDtree_filter(data, kdtree2d, target, truncated_radius,"polar")
        Sigmai = wrap_surfdens(kdtf_data2d,target)
        buffer_dict = wrap_quant2d(kdtf_data2d,target, Sigmai)
        for column_name in column_names
            Result_storage_dict[column_name][i] = buffer_dict[column_name]
        end
        Result_storage_dict["∇Sigmax"][i],Result_storage_dict["∇Sigmay"][i] = wrap_grad_surfdens(kdtf_data2d,target)
        Result_storage_dict["Sigma"][i] = Sigmai
    end
    @info "End 2d intepolate analysis."
    return Result_grid_dict
end

export density_analysis, surface_density_analysis, midH_analysis, pitch_analysis, pitch_other_analysis