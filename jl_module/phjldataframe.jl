using DataFrames
using LinearAlgebra
using Statistics

#Classes setting
struct phjlRawDataFrame
    dfdata :: DataFrame
    params :: Dict
end

#Method
function get_dim(data::phjlRawDataFrame)
    return hasproperty(data.dfdata,"z") ? 3 : 1
end

function pickup_position(data::phjlRawDataFrame,particle_index::Int)
    position = Vector{Float64}(undef, 3)
    variable = ["x","y","z"]
    for (i,var) in enumerate(variable)
        position[i] = data.dfdata[particle_index,var]
    end
    return position
end

function pickup_velocity(data::phjlRawDataFrame,particle_index::Int)
    velocity = Vector{Float64}(undef, 3)
    variable = ["vx","vy","vz"]
    for (i,var) in enumerate(variable)
        velocity[i] = data.dfdata[particle_index,var]
    end
    return velocity
end

function pickup_general_coordinate(data::phjlRawDataFrame,particle_index::Int)
    coordinate = Vector{Float64}(undef, 6)
    variable = String["x","y","z","vx","vy","vz"]
    for (i,var) in enumerate(variable)
        coordinate[i] = data.dfdata[particle_index,var]
    end
    return coordinate
end

function COM2primary(data_list, sinks_data:: phjlRawDataFrame)
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
    general_coordinateQ1 = pickup_general_coordinate(sinks_data,1)
    variable = ["x","y","z","vx","vy","vz"]
    for data in data_list
        for (i,var) in enumerate(variable)
            data.dfdata[:,var] .-= general_coordinateQ1[i]
        end
        data.params["COM_coordinate"] = -general_coordinateQ1
    end
end

# function primary2COM(data_list, sinks_data:: phjlRawDataFrame)
#     """
#     Transfer the coordinate into COM coordinate.
#     """
#     if (isa(data_list,Array))
#         nothing
#     elseif (isa(data_list,phjlRawDataFrame))
#         data_list = [data_list]
#     else
#         error("LoadError: Invaild Input in COM2primary")
#     end
#     if haskey(sinks_data.params,"COM_coordinate")
#         general_coordinateQ1 = data.params["COM_coordinate"]

#     end

#     general_coordinateQ1 = pickup_general_coordinate(sinks_data,1)
#     variable = ["x","y","z","vx","vy","vz"]
#     for data in data_list

#         for (i,var) in enumerate(variable)
#             data.dfdata[:,var] .-= general_coordinateQ1[i]
#         end
#     end
# end

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

function add_angular_momentum(data::phjlRawDataFrame)
    """add the angluar momentum w.r.t the current origin"""
    data.dfdata[!,"lx"] = (data.dfdata[!,"y"] .* data.dfdata[!,"vz"]) .- (data.dfdata[!,"z"] .* data.dfdata[!,"vy"])
    data.dfdata[!,"ly"] = (data.dfdata[!,"z"] .* data.dfdata[!,"vx"]) .- (data.dfdata[!,"x"] .* data.dfdata[!,"vz"])
    data.dfdata[!,"lz"] = (data.dfdata[!,"x"] .* data.dfdata[!,"vy"]) .- (data.dfdata[!,"y"] .* data.dfdata[!,"vx"])
end

function add_cylindrical(data::phjlRawDataFrame)
    if !(hasproperty(data.dfdata,"rnorm")) || !(hasproperty(data.dfdata,"vnorm"))
        add_norm(data)
        rename!(data.dfdata,"rnorm" => "s")
    end
    data.dfdata[!,"phi"] = atan.(data.dfdata[!,"y"],data.dfdata[!,"x"])
    sinphi = sin.(data.dfdata[!,"phi"])
    cosphi = cos.(data.dfdata[!,"phi"])
    data.dfdata[!,"vs"] = (cosphi.*data.dfdata[!,"vx"] + sinphi.*data.dfdata[!,"vy"])
    data.dfdata[!,"vphi"] = (cosphi.*data.dfdata[!,"vy"] - sinphi.*data.dfdata[!,"vx"])
end
