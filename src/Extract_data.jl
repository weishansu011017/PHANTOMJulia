abstract type Abstract_analysis_buffer end
abstract type Abstract_analysis end

struct Pitch_analysis_buffer <: Abstract_analysis_buffer
    time :: Float64
    data_dict :: Dict{Int, gridbackend}
    theta :: LinRange
    radius :: LinRange
    column_names :: Dict{Int,String}
end

mutable struct Pitch_analysis <: Abstract_analysis
    time :: Float64
    data_dict :: Dict{Int, Array{Float64}}
    theta :: Vector{Float64}
    radius :: Vector{Float64}
    column_names :: Dict{Int,String}
    params :: Dict{String,Any}
end

struct Spiral_analysis_buffer <: Abstract_analysis_buffer
    time :: Float64
    radius :: Float64
    data_dict :: Dict{Int, gridbackend}
    theta :: LinRange
    column_names :: Dict{Int,String}
end

mutable struct Spiral_analysis <: Abstract_analysis
    time :: Float64
    radius :: Float64
    data_dict :: Dict{Int, Array{Float64}}
    theta :: Vector{Float64}
    column_names :: Dict{Int,String}
    params :: Dict{String,Any}
end

function extract_number(str::String)
    m = match(r"(\d+)", str)
    return m !== nothing ? m.match : ""
end

function convert_field(value)
    if typeof(value) <: Dict{Int, gridbackend}
        converted_value = Dict{Int, Array}()
        for key in keys(value)
            converted_value[key] = value[key].grid
        end
        return converted_value
    elseif typeof(value) <: LinRange
        return collect(value)
    else
        return value
    end
end

function create_column_names(keys_order::Vector{String}, suffixes::Vector)
    column_names = ["theta"]
    for key in keys_order
        for suffix in suffixes
            push!(column_names, string(key, "_", suffix))
        end
    end
    return column_names
end
function create_column_dict(column_names::Array{String})
    column_format = Dict{Int, String}()

    total_length = 16
    index_length = 2

    for (i, name) in enumerate(column_names)
        index_str = lpad(i, index_length, '0')
        space_padding = total_length - length(index_str) - length(name) - 3 
        formatted_name = "[" * index_str * " " * repeat(" ", space_padding) * name * "]"
        column_format[i] = formatted_name
    end
    return column_format
end

function Pitch_analysis_buffer(time::Float64, Gas_grid::Dict, Dust_grid::Dict, key_values::Vector{String})
    column_names = create_column_names(key_values,["g","d"])
    column_dict = create_column_dict(column_names)
    gridbackend_dict = Dict{Int, gridbackend}()

    for key in key_values
        keyg = string(key, "_", "g")
        keyd = string(key, "_", "d")
        for dictkey in keys(column_dict)
            if occursin(keyg,column_dict[dictkey])
                try
                    gridbackend_dict[dictkey] = Gas_grid[key]
                catch e 
                    if isa(e,KeyError)
                        error("ExtractError: Missing key $(key) in Gas_grid!")
                    else
                        error("ExtractError: UnknownError")
                    end
                end
            end
            if occursin(keyd,column_dict[dictkey])
                try
                    gridbackend_dict[dictkey] = Dust_grid[key]
                catch e 
                    if isa(e,KeyError)
                        error("ExtractError: Missing key $(key) in Dust_grid!")
                    else
                        error("ExtractError: Error ($e) occured when loading Dust_grid!")
                    end
                end
            end
        end
    end

    radius = Gas_grid["Sigma"].axis[1]
    theta = Gas_grid["Sigma"].axis[2]
    return Pitch_analysis_buffer(
    time,
    gridbackend_dict,
    theta,
    radius,
    column_dict
    )
end

function Spiral_analysis_buffer(time::Float64,radius::Float64, Gas_grid::Dict, Dust_grid::Dict,key_values::Vector{String})
    column_names = create_column_names(key_values,["g","d"])
    column_dict = create_column_dict(column_names)
    gridbackend_dict = Dict{Int, gridbackend}()

    for key in key_values
        keyg = string(key, "_", "g")
        keyd = string(key, "_", "d")
        for dictkey in keys(column_dict)
            if occursin(keyg,column_dict[dictkey])
                try 
                    gridbackend_dict[dictkey] = Gas_grid[key]
                catch e
                    if isa(e,KeyError)
                        error("ExtractError: Missing key $(key) in Gas_grid!")
                    else
                        error("ExtractError: Error ($e) occured when loading Dust_grid!")
                    end
                end
            end
            if occursin(keyd,column_dict[dictkey])
                try 
                    gridbackend_dict[dictkey] = Dust_grid[key]
                catch
                    if isa(e,KeyError)
                        error("ExtractError: Missing key $(key) in Dust_grid!")
                    else
                        error("ExtractError: Error ($e) occured when loading Dust_grid!")
                    end
                end
            end
        end
    end

    theta = Gas_grid["Sigma"].axis[1]
    return Spiral_analysis_buffer(
    time,
    radius,
    gridbackend_dict,
    theta,
    column_dict
    )
end


function self_check_pitch(input :: Pitch_analysis_buffer)
    input_data = input.data_dict
    iaxis = Array{LinRange}(undef,2)
    iaxis[1] = input.radius
    iaxis[2] = input.theta
    for key in keys(input.column_names)
        column = input.column_names[key]
        if occursin("theta",column)
            continue
        end
        if typeof(input_data[key])<:gridbackend
            axis = input_data[key].axis
            for i in eachindex(iaxis)
                if !(iaxis[i] == axis[i])
                    error("AxisError: Mismatching of axis in $column")
                end
            end
        end
    end
    @info "Axis self-checking of Pitch_analysis success! "
end
function self_check_spiral(input :: Spiral_analysis_buffer)
    input_data = input.data_dict
    iaxis = Array{LinRange}(undef,1)
    iaxis[1] = input.theta
    for key in keys(input.column_names)
        column = input.column_names[key]
        if occursin("theta",column)
            continue
        end
        if typeof(input_data[key])<:gridbackend
            axis = input_data[key].axis
            for i in eachindex(iaxis)
                if !(iaxis[i] == axis[i])
                    error("AxisError: Mismatching of axis in $column")
                end
            end
        end
    end
    @info "Axis self-checking of Spiral_analysis success! "
end

function buffer2output(buffer_struct::Abstract_analysis_buffer)
    buffer_type = typeof(buffer_struct)
    if (buffer_type <: Pitch_analysis_buffer)
        self_check_pitch(buffer_struct)
    elseif (buffer_type <: Spiral_analysis_buffer)
        self_check_spiral(buffer_struct)
    else
        @warn "The type of struct $buffer_type is not recorded! It may cause some problem."
    end
    output_type = Symbol(replace(string(buffer_type), "_buffer" => ""))

    fields = fieldnames(buffer_type)
    converted_values = [convert_field(getfield(buffer_struct, f)) for f in fields]

    return eval(output_type)(converted_values...,Dict{String,Bool}())
end

function Write_pitch_dat(filename::String, data::Pitch_analysis)
    @info "Starting writting the output file for pitch analysis."
    numberfile = extract_number(filename)
    output_filename = "pitch_phjl_" * numberfile * ".dat"
    open(output_filename, "w") do f
        for i in eachindex(data.radius)
            @printf(f,"# Analysis data at t = %15.10E \n",data.time)
            @printf(f,"# theta_rad = %15.10E \n",data.radius[i])
            sorted_values = [data.column_names[key] for key in sort(collect(keys(data.column_names)))]
            println(f,"# ", join(sorted_values, "   "))
            for j in eachindex(data.theta)
                @printf(f,"%18.10E ",data.theta[j])
                for column_index in 1:length(data.column_names)
                    if column_index == 1
                        continue
                    end
                    @printf(f,"%18.10E ",data.data_dict[column_index][i,j])
                end
                println(f,"")
            end
        end
    end
    @info "Finished writting the output file for pitch analysis."
end

function Write_spiral_dat(filename::String, data::Spiral_analysis)
    @info "Starting writting the output file for spiral analysis."
    numberfile = extract_number(filename)
    output_filename = "spiral_phjl_" * numberfile * ".dat"
    open(output_filename, "w") do f
        @printf(f,"# Analysis data at t = %15.10E \n",data.time)
        @printf(f,"# theta_rad = %15.10E \n",data.radius)
        sorted_values = [data.column_names[key] for key in sort(collect(keys(data.column_names)))]
        println(f,"# ", join(sorted_values, "   "))
        for j in eachindex(data.theta)
            @printf(f,"%18.10E ",data.theta[j])
            for column_index in 1:length(data.column_names)
                if column_index == 1
                    continue
                end
                @printf(f,"%18.10E ",data.data_dict[column_index][j])
            end
            println(f,"")
        end
    end
    @info "Finished writting the output file for spiral analysis."
end

function Write_H5DF(filename::String, data::Abstract_analysis)
    type = typeof(data)
    numberfile = extract_number(filename)
    if type <: Pitch_analysis
        output_filename = "pitch_phjl_" * numberfile * ".h5"
    elseif type <: Spiral_analysis
        output_filename = "spiral_phjl_" * numberfile * ".h5"
    else
        error("TypeError: Non-supported type of data.")
    end
    type_data = string(type)
    h5open(output_filename,"w") do f
        write(f, "struct_type", type_data)
        for name in fieldnames(type)
            if String(name) != "params"
                val = getfield(data,name)
                if typeof(val) <: AbstractArray
                    write(f, string(name), val)
                elseif typeof(val) <: Dict
                    g = create_group(f, string(name))
                    for (key, value) in val
                        write(g, string(key), value)
                    end
                else
                    write(f, string(name), val)
                end
            end
        end
    end
end

function Write_output(filename::String, data::Abstract_analysis,mode::String, output_datatype::String="H5DF")
    if output_datatype=="H5DF"
        Write_H5DF(filename,data)
    elseif output_datatype=="dat"
        if mode=="Spiral"
            Write_spiral_dat(filename,data)
        elseif mode=="Pitch"
            Write_pitch_dat(filename,data)
        else
            @warn "OutputWarn: The output mode $mode is unsupported. Changing the output datatype to H5DF!"
            Write_H5DF(filename,data)
        end
    else
        error("OutputError: Unsupported output datatype $(output_datatype)!")
    end
end

for name in names(@__MODULE__; all=true)
    if name ∉ (:include, :eval) && isdefined(@__MODULE__, name)
        @eval export $name
    end
end


