function add_necessary_quantity(data::phjlRawDataFrame)
    if !(hasproperty(data.dfdata,"rho"))
        add_rho(data)
    end
end

function _standard_estimate_h_intepolate(dfdata::DataFrame, rnorm::Vector,  distance_threshold::Float64)
    """Assume data has been transfered to the reference_point-based coordinate"""
    indices = findall(x -> x <= 1.2*distance_threshold, rnorm)
    filtered_df = dfdata[indices,:]
    fallback :: Float32 = 0.0
    h_intepolate :: Float32 = 0.0
    total_weight :: Float32 = 0.0
    j=0
    for particle in eachrow(filtered_df)
        j += 1
        distance = rnorm[j]
        h_particle = particle.h
        if distance < distance_threshold
            return h_particle
        else
            weight = 1.0/distance
            h_intepolate += weight * h_particle
            total_weight += weight
        end
    end
    return total_weight>0 ? h_intepolate/total_weight : fallback

end


function _easy_estimate_h_intepolate(dfdata::DataFrame, rnorm::Vector, distance_threshold::Float64)
    """
    Give a specific smoothed radius to calculate
    Assume data has been transfered to the reference_point-based coordinate
    """
    indices = findall(x -> x <= 1.2*distance_threshold, rnorm)
    if (isempty(indices))
        return Float32(0.0)
    end
    min_index = indices[argmin(rnorm[indices])]
    return dfdata[min_index,"h"]
end

function estimate_h_intepolate(data::phjlRawDataFrame, rnorm::Vector, distance_threshold::Float64,mode::String ="intep")
    """
    Give a specific smoothed radius to calculate
    Assume data has been transfered to the reference_point-based coordinate

    mode: 
    "intep": Intepolate h by standard.
    "closest": Choose h of the closest particles.
    "mean": use the mean value of h
    """
    dfdata = data.dfdata
    if !(haskey(data.params,"h_mean"))
        add_mean_h(data)
    end
    if (mode == "mean")
        return data.params["h_mean"]
    elseif (mode == "closest")
        return _easy_estimate_h_intepolate(dfdata, rnorm, distance_threshold)
    elseif (mode == "intep")
        return _standard_estimate_h_intepolate(dfdata, rnorm, distance_threshold)
    else
        error("IntepolateError: Invaild mode of h calculattion.")
    end
end

function density(data::phjlRawDataFrame, reference_point::Vector, smoothed_kernal:: Function = M4_spline,h_mode::String="intep", kind_flag::String = "cart")
    """
    Here recommended to use a single type of particle.
    kind_flag is the coordinate system that the reference_point is given
    reference_point is in "3D"
    "cart" = cartitian
    "polar" = cylindrical
    """
    if kind_flag == "polar"
        reference_point = _cylin2cart(reference_point)
    end
    truncate_multiplier = KernelFunctionValid()[nameof(smoothed_kernal)]
    rnorm = get_rnorm_ref(data, reference_point)
    particle_mass = data.params["mass"]
    h_intepolate = estimate_h_intepolate(data,rnorm,truncate_multiplier,h_mode) #_easy_estimate_h_intepolate(dfdata, rnorm, 1.0)
    if h_intepolate == 0.0
        density = 0.0
    else
        truncate_radius = truncate_multiplier * h_intepolate
        filtered_rnorm = filter(r -> r <= truncate_radius, rnorm)
        density = sum(particle_mass.*Smoothed_kernel_function.(smoothed_kernal,h_intepolate,filtered_rnorm,3))
    end
    return density
end

function surface_density(data::phjlRawDataFrame, reference_point::Vector, smoothed_kernal:: Function = M4_spline,h_mode::String="intep", kind_flag::String = "cart")
    """
    Here recommended to use a single type of particle.
    kind_flag is the coordinate system that the reference_point is given
    reference_point is in "2D"
    "cart" = cartitian
    "polar" = polar
    """
    if kind_flag == "polar"
        reference_point = _cylin2cart(reference_point)
    end
    truncate_multiplier = KernelFunctionValid()[nameof(smoothed_kernal)]
    snorm = get_snorm_ref(data, reference_point)
    particle_mass = data.params["mass"]
    h_intepolate = estimate_h_intepolate(data,snorm,truncate_multiplier,h_mode) #_easy_estimate_h_intepolate(dfdata, rnorm, 1.0)
    if h_intepolate == 0.0
        surface_density = 0.0
    else
        truncate_radius = truncate_multiplier * h_intepolate
        filtered_snorm = filter(r -> r <= truncate_radius, snorm)
        surface_density = sum(particle_mass.*Smoothed_kernel_function.(smoothed_kernal,h_intepolate,filtered_snorm,2))
    end
    return surface_density
end

function gradient_surface_density(data::phjlRawDataFrame, reference_point::Vector, smoothed_kernal:: Function = M4_spline,h_mode::String="intep", kind_flag::String = "cart")
    """
    Here recommended to use a single type of particle.
    kind_flag is the coordinate system that the reference_point is given
    reference_point is in "2D"
    "cart" = cartitian
    "polar" = polar
    """
    if kind_flag == "polar"
        reference_point = _cylin2cart(reference_point)
    end
    truncate_multiplier = KernelFunctionValid()[nameof(smoothed_kernal)]
    snorm,xyref = get_s_ref(data, reference_point)
    particle_mass = data.params["mass"]
    h_intepolate = estimate_h_intepolate(data,snorm,truncate_multiplier,h_mode) #_easy_estimate_h_intepolate(dfdata, rnorm, 1.0)
    grad_surface_density = zeros(Float64,2)
    if h_intepolate == 0.0
        return grad_surface_density
    else
        truncate_radius = truncate_multiplier * h_intepolate
        mask_snorm = snorm .< truncate_radius
        xy_filtered = xyref[mask_snorm,:]
        buffer_array = zeros(Float64,size(xy_filtered))
        for i in 1:size(xy_filtered)[1]
            buffer_array[i,:] = particle_mass.*Smoothed_greident_kernel_function(smoothed_kernal,h_intepolate,xy_filtered[i,:])
        end
        for j in eachindex(grad_surface_density)
            grad_surface_density[j] = sum(buffer_array[:,j])
        end
        return grad_surface_density
    end
end

function quantity_intepolate_2D(data::phjlRawDataFrame, reference_point::Vector,Sigmai::Float64, column_names::Vector{String}, smoothed_kernal:: Function = M4_spline,h_mode::String="intep", kind_flag::String = "cart")
    """
    Here recommended to use a single type of particle.
    kind_flag is the coordinate system that the reference_point is given
    reference_point is in "2D"
    "cart" = cartitian
    "polar" = polar
    """
    for column_name in column_names
        if !(hasproperty(data.dfdata,column_name))
            error("IntepolateError: No matching column '$(column_name)'.")
        end
    end

    if kind_flag == "polar"
        reference_point = _cylin2cart(reference_point)
    end
    
    quantity_result = Dict{String,Float64}()
    truncate_multiplier = KernelFunctionValid()[nameof(smoothed_kernal)]
    snorm = get_snorm_ref(data, reference_point)
    particle_mass = data.params["mass"]
    h_intepolate = estimate_h_intepolate(data, snorm, truncate_multiplier,h_mode)
    dfdata = data.dfdata
    

    if h_intepolate == 0.0
        density = 0.0
    else
        truncate_radius = truncate_multiplier * h_intepolate
        indices = findall(x -> x <= truncate_radius, snorm)
        filtered_dfdata = dfdata[indices, :]
        filtered_snorm = filter(r -> r <= truncate_radius, snorm)
        for column_name in column_names
            filtered_dfdata[!,column_name] ./= Sigmai
            quantity_result[column_name] = sum(particle_mass.*(filtered_dfdata[!,column_name]).*(Smoothed_kernel_function.(smoothed_kernal,h_intepolate,filtered_snorm,2)))
        end
    end
    return quantity_result
end

function quantity_intepolate(data::phjlRawDataFrame, reference_point::Vector, column_names::Vector{String}, smoothed_kernal:: Function = M4_spline,h_mode::String="intep", kind_flag::String = "cart")
    """
    Here recommended to use a single type of particle.
    kind_flag is the coordinate system that the reference_point is given
    "cart" = cartitian
    "polar" = cylindrical
    """
    for column_name in column_names
        if !(hasproperty(data.dfdata,column_name))
            error("IntepolateError: No matching column '$(column_name)'.")
        end
    end
    if !(hasproperty(data.dfdata,"rho"))
        add_rho(data)
    end
    if kind_flag == "polar"
        reference_point = _cylin2cart(reference_point)
    end
    quantity_result = Dict{String,Float64}()
    truncate_multiplier = KernelFunctionValid()[nameof(smoothed_kernal)]
    rnorm = get_rnorm_ref(data, reference_point)
    particle_mass = data.params["mass"]
    h_intepolate = estimate_h_intepolate(data, rnorm, truncate_multiplier,h_mode)
    dfdata = data.dfdata
    

    if h_intepolate == 0.0
        density = 0.0
    else
        truncate_radius = truncate_multiplier * h_intepolate
        indices = findall(x -> x <= truncate_radius, rnorm)
        filtered_dfdata = dfdata[indices, :]
        filtered_rnorm = filter(r -> r <= truncate_radius, rnorm)
        for column_name in column_names
            filtered_dfdata[!,column_name] ./= filtered_dfdata[!,"rho"]
            quantity_result[column_name] = sum(particle_mass.*(filtered_dfdata[!,column_name]).*(Smoothed_kernel_function.(smoothed_kernal,h_intepolate,filtered_rnorm,3)))
        end
    end
    return quantity_result
end

