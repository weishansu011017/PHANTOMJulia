include("mathematical_tools.jl")
include("kernel_function.jl")
include("phjldataframe.jl")
include("grid.jl")



function _estimate_h_intepolate(dfdata::DataFrame, h_threshold::Float64, distance_threshold::Float64)
    """Assume data has been transfered to the reference_point-based coordinate"""
    filtered_df = dfdata[dfdata.rnorm .< distance_threshold,:]

    fallback :: Float32 = 0.0
    h_intepolate :: Float32 = 0.0
    total_weight :: Float32 = 0.0

    for particle in eachrow(filtered_df)
        distance = particle.rnorm
        h_particle = particle.h
        if distance < h_threshold
            return h_particle
        else
            weight = 1.0/distance
            h_intepolate += weight * h_particle
            total_weight += weight
        end
    end
    return total_weight>0 ? h_intepolate/total_weight : fallback

end

function _easy_estimate_h_intepolate(dfdata::DataFrame, h_threshold::Float64, distance_threshold::Float64)
    """Assume data has been transfered to the reference_point-based coordinate"""
    min_index = argmin(dfdata[!,"rnorm"])
    if dfdata[min_index,"rnorm"] > distance_threshold
        return Float32(0.0)
    else
        return dfdata[min_index,"h"]
    end
end

function density(data::phjlRawDataFrame, reference_point::Vector, smoothed_kernal:: Function = M4_spline, kind_flag::String = "cart")
    """
    Here recommended to use a single type of particle.
    kind_flag is the coordinate system that the reference_point is given
    "cart" = cartitian
    "cylin" = cylindrical
    """
    if kind_flag == "cylin"
        reference_point = cylin2cart(reference_point)
    end
    if !(hasproperty(data.dfdata,"rho"))
        add_rho(data)
    end
    dfdata = copy(data.dfdata)

    for (i,dir) in enumerate(String["x","y","z"])
        dfdata[:,dir] .-= reference_point[i]  
    end
    add_norm(dfdata)
    particle_mass = data.params["massoftype"]
    h_intepolate = _estimate_h_intepolate(dfdata, 0.1, 1.0)
    if h_intepolate == 0.0
        density = 0.0
    else
        density = sum(particle_mass.*Smoothed_kernel_function.(smoothed_kernal,h_intepolate,dfdata[!,"rnorm"]))
    end
    return density
end


function quantity_intepolate(data::phjlRawDataFrame, reference_point::Vector, target::Vector, smoothed_kernal:: Function = M4_spline, kind_flag::String = "cart")
    """
    Here recommended to use a single type of particle.
    kind_flag is the coordinate system that the reference_point is given
    "cart" = cartitian
    "cylin" = cylindrical

    target: target columns that want to be intepolate
    ["rnorm"] or ["lx","ly","lz"]
    """
    if kind_flag == "cylin"
        reference_point = cylin2cart(reference_point)
    end
    if !(hasproperty(data.dfdata,"rho"))
        add_rho(data)
    end
    for column in target
        if !(hasproperty(data.dfdata,column))
            error("IntepolatingError: Missing column properties")
        end
    end

    datacopy = deepcopy(data)
    for (i,dir) in enumerate(String["x","y","z"])
        datacopy.dfdata[:,dir] .-= reference_point[i]  
    end

    particle_mass = datacopy.params["massoftype"]
    h_intepolate = _easy_estimate_h_intepolate(datacopy, 0.1, 2.0)
    result = zeros(length(target))
    if h_intepolate != 0.0
        for (i,column) in enumerate(target)
            result[i] = sum(particle_mass.*(datacopy.dfdata[!,column]./datacopy.dfdata[!,"rho"]).*Smoothed_kernel_function.(smoothed_kernal,h_intepolate,datacopy.dfdata[!,"rnorm"]))
        end
    end
    return result
end

function surface_density(data::phjlRawDataFrame, reference_point::Vector, Hovr::Float64, smoothed_kernal:: Function = M4_spline, z_seperate::Int = 4,kind_flag::String = "cylin")
    """
    Calculate the surface density at (x,y) or (r,θ) in the inteval z ∈ [Hmin,Hmax] 
    """
    function wrap_density(point)
        return density(data,point,smoothed_kernal,kind_flag)
    end
    if kind_flag == "cart"
        reference_point = cart2cylin(reference_point)
    end
    r = reference_point[1]
    theta = reference_point[2]
    println("Dealing with $r and $theta")
    Hlim = r * Hovr
    zrange = LinRange(-Hlim, Hlim, z_seperate)
    point_array = [[r,theta,z] for z in zrange]
    rho = zeros(Float64,z_seperate)
    for i in eachindex(point_array)
        rho[i] = wrap_density(point_array[i])
    end
    rho = wrap_density.(point_array)
    surface_density,error = Integral_1d(zrange,rho,[-Hlim,Hlim])
    return surface_density
end

