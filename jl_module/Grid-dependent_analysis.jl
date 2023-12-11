include("grid.jl")
include("phjldataframe.jl")




function surface_density_analysis(data::phjlRawDataFrame,r_params::Tuple{Float64,Float64,Int}, theta_n :: Int, Hovr::Float64, smoothed_kernal:: Function = M4_spline, z_seperate::Int = 20)
    """
    Calculate the surface density in a whole disc
    r_params = (rmin,rmax,rn)
    """
    function wrap_surd(point)
        return surface_density(data, point, Hovr, smoothed_kernal,z_seperate)
    end
    imin = [r_params[1], 0.0]
    imax = [r_params[2], 2*pi]
    in = [r_params[3],theta_n]
    disc_grid = disc_2d_grid_generator(imin,imax,in)
    gridv = generate_coordinate_grid(disc_grid)
    gridv_line = gridv[:]
    storage = zeros(Float64, in[1]*in[2])
    # gridv = broadcast(point -> surface_density(data, point, Hovr, smoothed_kernal,z_seperate), gridv)
    @threads for i in eachindex(gridv_line)
        storage[i] = wrap_surd(gridv_line[i])
    end
    gridv = reshape(storage, in...)
    return gridv
end
