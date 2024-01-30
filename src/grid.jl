struct gridbackend
    grid :: Array
    axis :: Vector{LinRange}
    dimension :: Vector{Int}
end

function meshgrid(arrays::AbstractVector...)
    nd = length(arrays)
    grids = [repeat(reshape(a, ntuple(d -> d == i ? length(a) : 1, nd)...), 
                    ntuple(d -> d == i ? 1 : length(arrays[d]), nd)...) for (i, a) in enumerate(arrays)]
    return grids
end

function meshgrid(T::gridbackend)
    iaxis = T.axis
    meshgrids = meshgrid(iaxis...)
    return meshgrids
end

function generate_empty_grid(imin::Vector{Float64}, imax::Vector{Float64}, dimension::Vector{Int}, type::Type = Float64)
    """
    Generate a empty grid.
    imin: minimum of each axis e.g. [xmin,ymin,zmin]
    imax: maximum of each axis e.g. [xmax,ymax,zmax]
    dimension: The dimension of the grid (how much do the axies been seperated.) e.g 3x3x4 matrix => [3,3,4]
    type:: Type in the data e.g. Float64

    return grid:: gridbackend
    """
    if size(imin) == size(imax) == size(dimension)
        nothing
    else
        error("GridGeneratingError: Illegal input value.")
    end
    iaxis = Array{LinRange}(undef,length(dimension))
    for (i,num) in enumerate(dimension)
        iaxis[i] = LinRange(imin[i],imax[i],num)
    end
    grid :: Array = zeros(type, dimension...)
    return gridbackend(grid,iaxis,dimension)
end

function coordinate(T :: gridbackend, element :: Tuple)
    """
    Returning the corresponding coordinate of specific element.
    """
    if length(T.dimension) != length(element)
        error("GridLodingError: Mismatching of dimension between input element and grid.")
    end

    result :: Array = zeros(Float64,length(T.dimension))
    for i in eachindex(result)
        result[i] = T.axis[i][element[i]]
    end
    return result 
end

function generate_coordinate_grid(T :: gridbackend)
    """
    Returning an Array which contain the coordinate of the coorespoing element.
    """
    dims = T.dimension

    coordinates_array = Array{Array{Float64,1}, length(dims)}(undef, dims...)
    for idx in CartesianIndices(coordinates_array)
        element_index = Tuple(idx)
        coordinates_array[idx] = coordinate(T, element_index)
    end
    return coordinates_array
end

function disc_2d_grid_generator(imin::Vector ,imax::Vector, in::Vector)
    """
    Generate a disk grid 
    imin = [rmin,θmin]
    imax = [rmax,θmax]
    in = [rn,θn]

    return: grid, iaxis
    grid :: Array = zeros(rn,θn)
    iaxis :: Vector = [
        for (i,num) in enumerate(in)
            LinRange(...)
        end
    ] 
    """
    if size(imin) == size(imax) == size(in)
        nothing
    else
        error("GridGeneratingError: Illegal input value.")
    end
    if (imax[2]>2*pi) || (imin[2]<0)
        error("GridGeneratingError: Illegal theta value.")
    end
    backend = generate_empty_grid(imin,imax,in)
    return backend
end
