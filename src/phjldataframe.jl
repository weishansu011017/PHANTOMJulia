#struct setting
struct phjlRawDataFrame
    dfdata :: DataFrame
    params :: Dict
end

"""
All of the function that can modify the phjlRawDataFrame. 
Can add any quantity for each particles for SPH intepolation.
"""

#Method
function get_dim(data::phjlRawDataFrame)
    return hasproperty(data.dfdata,"z") ? 3 : 2
end

function add_mean_h(data::phjlRawDataFrame)
    data.params["h_mean"] = mean(data.dfdata[!,"h"])
end

function get_truncated_radius(data::phjlRawDataFrame, poffset::Float64=0.5, smoothed_kernal:: Function = M4_spline)
    if !(haskey(data.params, "h_mean"))
        add_mean_h(data)
    end
    radius = (KernelFunctionValid()[nameof(smoothed_kernal)] + poffset)*data.params["h_mean"]
    return radius
end

function print_params(data::phjlRawDataFrame, pause::Bool=false)
    allkeys = sort(collect(keys(data.params)))
    for key in allkeys
        println("$(key) => $(data.params[key])")
    end
    if pause
        readline()
    end
end

function pickup_general_coordinate(data::phjlRawDataFrame,particle_index::Int)
    coordinate = Vector{Float64}(undef, 6)
    variable = String["x","y","z","vx","vy","vz"]
    for (i,var) in enumerate(variable)
        coordinate[i] = data.dfdata[particle_index,var]
    end
    return coordinate
end

function Generate_KDtree(data::phjlRawDataFrame,dim::Int)
    if (dim==2)
        position_array = hcat(data.dfdata[!,:x],data.dfdata[!,:y])'
    elseif (dim==3)
        position_array = hcat(data.dfdata[!,:x],data.dfdata[!,:y],data.dfdata[!,:z])'
    end
    kdtree = KDTree(position_array) 
    return kdtree
end

function KDtree_filter(data::phjlRawDataFrame, kdtree::KDTree, target::Vector, radius::Float64, kind_flag::String = "cart")
    """
    Here recommended to use a single type of particle.
    kind_flag is the coordinate system that the reference_point is given
    reference_point is in "2D"
    "cart" = cartitian
    "polar" = polar
    """
    if kind_flag == "polar"
        target = _cylin2cart(target)
    end
    dim = length(first(kdtree.data))
    if (dim != length(target))
        error("DimensionalError: The kdtree is constructed in $(dim)-d, but the given target is in $(length(target))-d.")
    end
    kdtf_dfdata = data.dfdata[inrange(kdtree, target, radius),:]
    kdtf_data = phjlRawDataFrame(kdtf_dfdata, data.params)
    return kdtf_data
end

function get_position_array(data::phjlRawDataFrame)
    position_array = hcat(data.dfdata[!, :x], data.dfdata[!, :y], data.dfdata[!, :z])
    return position_array
end

function get_rnorm_ref(data::phjlRawDataFrame, reference_position::Vector{Float64})
    xt,yt,zt = reference_position
    x, y, z = data.dfdata[!,"x"], data.dfdata[!,"y"], data.dfdata[!,"z"]
    rnorm = sqrt.((x .- xt).^2 + (y .- yt).^2 + (z.-zt).^2)
    return rnorm
end

function get_snorm_ref(data::phjlRawDataFrame, reference_position::Vector{Float64})
    if length(reference_position) == 2
        xt,yt = reference_position
    elseif length(reference_position) == 3
        xt,yt,zt = reference_position
    else
        error("DimensionalError: Wrong length for reference_position.")
    end
    x, y = data.dfdata[!,"x"], data.dfdata[!,"y"]
    snorm = sqrt.((x .- xt).^2 + (y .- yt).^2)
    return snorm
end

function get_snorm(data::phjlRawDataFrame)
    return get_snorm_ref(data,[0.0,0.0,0.0])
end

function get_s_ref(data::phjlRawDataFrame,reference_position::Vector{Float64})
    if length(reference_position) == 2
        xt,yt = reference_position
    elseif length(reference_position) == 3
        xt,yt,zt = reference_position
    else
        error("DimensionalError: Wrong length for reference_position.")
    end
    x = xt .- data.dfdata[!,"x"] 
    y = yt .- data.dfdata[!,"y"] 
    xy = hcat(x,y)
    snorm = sqrt.(x.^2 + y.^2)
    return snorm, xy

end

function COM2primary(data_list::Array, sinks_data:: phjlRawDataFrame,sink_particle_id::Int)
    """
    Transfer the coordinate into a primary star-based coordinate.
    """
    if (isa(data_list,Array))
        nothing
    elseif (isa(data_list,phjlRawDataFrame))
        data_list = [data_list]
    else
        error("LoadError: Invaild Input in COM2primary")
    end
    general_coordinateQ1 = pickup_general_coordinate(sinks_data,sink_particle_id)
    variable = ["x","y","z","vx","vy","vz"]
    for data in data_list
        for (i,var) in enumerate(variable)
            data.dfdata[:,var] .-= general_coordinateQ1[i]
        end
        data.params["COM_coordinate"] .-= general_coordinateQ1
        data.params["Origin_located"] = "primary"
        data.params["Origin_sink_mass"] = sinks_data.dfdata[sink_particle_id,"m"]
    end
end

function add_cylindrical(data::phjlRawDataFrame)
    data.dfdata[!,"s"] = sqrt.(data.dfdata[!,"x"].^2+data.dfdata[!,"y"].^2)
    data.dfdata[!,"theta"] = atan.(data.dfdata[!,"y"],data.dfdata[!,"x"])
    sintheta = sin.(data.dfdata[!,"theta"])
    costheta = cos.(data.dfdata[!,"theta"])
    data.dfdata[!,"vr"] = (costheta.*data.dfdata[!,"vx"] + sintheta.*data.dfdata[!,"vy"])
    data.dfdata[!,"vtheta"] = (costheta.*data.dfdata[!,"vy"] - sintheta.*data.dfdata[!,"vx"])
end

function add_norm(data::phjlRawDataFrame)
    data.dfdata[!,"vnorm"] = sqrt.(data.dfdata[!,"vx"].^2 + data.dfdata[!,"vy"].^2 + data.dfdata[!,"vz"].^2)
    data.dfdata[!,"rnorm"] = sqrt.(data.dfdata[!,"x"].^2 + data.dfdata[!,"y"].^2 + data.dfdata[!,"z"].^2)
end

function add_norm(dfdata::DataFrame)
    dfdata[!,"vnorm"] = sqrt.(dfdata[!,"vx"].^2 + dfdata[!,"vy"].^2 + dfdata[!,"vz"].^2)
    dfdata[!,"rnorm"] = sqrt.(dfdata[!,"x"].^2 + dfdata[!,"y"].^2 + dfdata[!,"z"].^2)
end

function add_rho(data::phjlRawDataFrame)
    particle_mass = data.params["massoftype"]
    hfact = data.params["hfact"]
    d = get_dim(data)
    data.dfdata[!,"rho"] = particle_mass .* (hfact./data.dfdata[!,"h"]).^(d) 
end

function add_Kepelarian_azimuthal_vecocity(data::phjlRawDataFrame)
    if !(hasproperty(data.dfdata,"s"))
        add_cylindrical(data)
    end
    G = 1.0
    M = data.params["Origin_sink_mass"]
    data.dfdata[!,"vtheta_k"] = sqrt.((G*M)./data.dfdata[!,"s"])
    data.dfdata[!,"vtheta_sub"] = data.dfdata[!,"vtheta"] - data.dfdata[!,"vtheta_k"]
    data.dfdata[!,"vsubnorm"] = sqrt.(data.dfdata[!,"vr"].^2 + data.dfdata[!,"vtheta_sub"].^2 + data.dfdata[!,"vz"].^2)
end

function add_kinetic_energy(data::phjlRawDataFrame)
    """
    Add the Kinetic energy in current frame.
    """
    if !(hasproperty(data.dfdata, "vnorm"))
        add_norm(data)
    end
    dfdata = data.dfdata
    particle_mass = data.params["massoftype"]
    data.dfdata[!,"KE"] = (particle_mass/2).*dfdata[!,"vnorm"]
end

function add_specialized_angular_momentum(data::phjlRawDataFrame)
    """add the angluar momentum w.r.t the current origin"""
    data.dfdata[!,"lx"] = (data.dfdata[!,"y"] .* data.dfdata[!,"vz"]) .- (data.dfdata[!,"z"] .* data.dfdata[!,"vy"])
    data.dfdata[!,"ly"] = (data.dfdata[!,"z"] .* data.dfdata[!,"vx"]) .- (data.dfdata[!,"x"] .* data.dfdata[!,"vz"])
    data.dfdata[!,"lz"] = (data.dfdata[!,"x"] .* data.dfdata[!,"vy"]) .- (data.dfdata[!,"y"] .* data.dfdata[!,"vx"])
    data.dfdata[!,"lnorm"] = sqrt.(data.dfdata[!,"lx"].^2 + data.dfdata[!,"ly"].^2 + data.dfdata[!,"lz"].^2)
end

function add_disc_normalized_angular_momentum(data::phjlRawDataFrame, rmin::Float64, rmax::Float64)
    """calculate the disc angular momentum"""
    if !(hasproperty(data.dfdata,"lx")) || !(hasproperty(data.dfdata,"ly")) || !(hasproperty(data.dfdata,"lz"))
        add_specialized_angular_momentum(data)
    end
    snorm = get_snorm(data)
    ldisc = zeros(Float64,3)
    disc_particles = (snorm.>rmin).&(snorm.<rmax)
    for (i,dir) in enumerate(["lx","ly","lz"])
        ldisc[i] = mean(data.dfdata[disc_particles,dir])
    end
    ldisc ./= norm(ldisc)
    data.params["ldisc"] = ldisc
end

function add_tilt(data::phjlRawDataFrame, rmin::Float64, rmax::Float64)
    if !(hasproperty(data.dfdata,"lx")) || !(hasproperty(data.dfdata,"ly")) || !(hasproperty(data.dfdata,"lz"))
        add_disc_normalized_angular_momentum(data,rmin,rmax)
    end
    if !(hasproperty(data.dfdata,"rnorm"))
        add_norm(data)
    end
    rlproject = (data.dfdata[!,"x"].*data.dfdata[!,"lx"] + data.dfdata[!,"y"].*data.dfdata[!,"ly"] + data.dfdata[!,"z"].*data.dfdata[!,"lz"])./data.dfdata[!,"lnorm"]
    nonzero_rnorm = data.dfdata[!,"rnorm"] .!= 0
    data.dfdata[!,"tilt"] = asin.(rlproject[nonzero_rnorm]./data.dfdata[nonzero_rnorm,"rnorm"])
end

function rotate_to_primary_L(data_list::Array, rmin::Float64, rmax::Float64, target_laxis::Union{Nothing, Vector{Float64}} = nothing)
    """
    Rotate the whole data to make z become angular_momentum_vector.
    Will take the angular momentum information from the first file. 

    if no (data_list[1].params['ldisc']) => Make it

            1      0      0
    Rx = [  0   cos(θx) -sin(θx)]
            0   sin(θx)  cos(θx) 

          cos(θy)  0     sin(θy)
    Ry = [  0      1      0     ]
         -sin(θy)  0     cos(θy)

    R = RxRy => rotate y axis and then x axis
    l = (lx,ly,lz), lxz = N(lx,0,lz), N = 1/√lx^2 + lz^2
    θy = Nlz, θx = N(lx^2 + lz^2) = 1/N
    """
    for data in data_list
        if (data.params["Origin_located"] == "primary")
            nothing
        elseif (data.params["Origin_located"] == "COM")
            COM2primary(data,data_list[end],1)
        else
            error("OriginLocatedError: Wrong origin located.")
        end
    end
    if isnothing(target_laxis)
        if !(haskey(data_list[1].params,"ldisc"))
            add_disc_normalized_angular_momentum(data_list[1],rmin,rmax)
            laxis = data_list[1].params["ldisc"]
        end
    else
        laxis = target_laxis
    end
    if laxis[3] < 0
        laxis = -laxis
    end
    lx,ly,lz = laxis
    N = 1/sqrt(lx^2 + lz^2)
    costhy = N*lz
    sinthy = sin(acos(costhy))
    costhx = 1/N
    sinthx = sin(acos(costhx))
    sinsinxy = sinthx * sinthy
    cossinxy = costhx * sinthy
    sincosxy = sinthx * costhy
    coscosxy = costhx * costhy
    for data in data_list
        dfdata = data.dfdata
        copydfdata = deepcopy(dfdata)
        dfdata[!, "x"] = costhy*copydfdata[!, "x"] + 0*copydfdata[!, "y"] + sinthy*copydfdata[!, "z"]
        dfdata[!, "y"] = sinsinxy*copydfdata[!, "x"] + costhx*copydfdata[!, "y"] - sincosxy*copydfdata[!, "z"]
        dfdata[!, "z"] = -cossinxy*copydfdata[!, "x"] + sinthx*copydfdata[!, "y"] + coscosxy*copydfdata[!, "z"]
        dfdata[!, "vx"] = costhy*copydfdata[!, "vx"] + 0*copydfdata[!, "vy"] + sinthy*copydfdata[!, "vz"]
        dfdata[!, "vy"] = sinsinxy*copydfdata[!, "vx"] + costhx*copydfdata[!, "vy"] - sincosxy*copydfdata[!, "vz"]
        dfdata[!, "vz"] = -cossinxy*copydfdata[!, "vx"] + sinthx*copydfdata[!, "vy"] + coscosxy*copydfdata[!, "vz"]
        add_disc_normalized_angular_momentum(data,rmin,rmax)
    end
    return laxis
end

function add_eccentricity(data::phjlRawDataFrame,sink_data::phjlRawDataFrame,G::Float64=1.0)
    """
    Add the eccentricity for each particle with respect to current origin. Default to be primary.
    """
    if !(haskey(sink_data.params, "Origin_located")) || (sink_data.params["Origin_located"] == "COM")
        error("OriginLocatedError: Wrong origin located.")
    elseif (sink_data.params["Origin_located"] == "primary") || (sink_data.params["Origin_located"] == "secondary")
        nothing
    else
        error("OriginLocatedError: The origin is not valid for the eccentricity calculation.")
    end
    M1 = sink_data.params["Origin_sink_mass"]
    μ = G * M1
    dfdata = data.dfdata
    if !(hasproperty(dfdata, "rnorm")) || !(hasproperty(dfdata, "vnorm"))
        add_norm(dfdata)
    end
    x, y, z = dfdata[!,"x"], dfdata[!,"y"], dfdata[!,"z"]
    vx, vy, vz = dfdata[!,"vx"], dfdata[!,"vy"], dfdata[!,"vz"]
    rnorm = dfdata[!,"rnorm"]
    vnorm = dfdata[!,"vnorm"]
    rdotv = (x.*vx) .+ (y.*vy) .+ (z.*vz)
    vnorm2 = vnorm .^2
    invrnorm = 1 ./rnorm
    ex = ((vnorm2./μ) .- invrnorm).* x - (rdotv./μ) .* vx
    ey = ((vnorm2./μ) .- invrnorm).* y - (rdotv./μ) .* vy
    ez = ((vnorm2./μ) .- invrnorm).* z - (rdotv./μ) .* vz
    dfdata[!,"e"] = sqrt.(ex.^2 + ey.^2 + ez.^2)
end

for name in names(PHANTOMJulia; all=true)
    if name ∉ (:include, :eval) && isdefined(PHANTOMJulia, name)
        @eval export $name
    end
end