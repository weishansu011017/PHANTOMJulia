struct Pitch_analysis_buffer
    time :: Float64
    rho_g :: gridbackend
    rho_d :: gridbackend
    vr_g :: gridbackend
    vr_d :: gridbackend
    vphi_g :: gridbackend
    vphi_d :: gridbackend
    sigma_gas :: gridbackend
    sigma_dust :: gridbackend
    theta :: LinRange
    radius :: LinRange
    column_names :: OrderedDict
end

struct Pitch_analysis
    time :: Float64
    rho_g :: Array
    rho_d :: Array
    vr_g :: Array
    vr_d :: Array
    vphi_g :: Array
    vphi_d :: Array
    sigma_gas :: Array
    sigma_dust :: Array
    theta :: Vector
    radius :: Vector
    column_names :: OrderedDict
end

struct Spiral_analysis_buffer
    time :: Float64
    radius :: Float64
    sigma_gas :: gridbackend
    tilt_gas :: gridbackend
    sigma_dust :: gridbackend
    tilt_dust :: gridbackend
    theta :: LinRange
    column_names :: OrderedDict
end

struct Spiral_analysis
    time :: Float64
    radius :: Float64
    sigma_gas :: Array
    tilt_gas :: Array
    sigma_dust :: Array
    tilt_dust :: Array
    theta :: Vector
    column_names :: OrderedDict
end

function extract_number(str::String)
    m = match(r"(\d+)", str)
    return m !== nothing ? m.match : ""
end

function convert_field(value)
    if typeof(value) <: gridbackend
        return value.grid  
    elseif typeof(value) <: LinRange
        return collect(value)
    else
        return value
    end
end

function create_column_dict(column_names::Array{String})
    column_format = OrderedDict{String, String}()

    total_length = 16
    index_length = 2

    for (i, name) in enumerate(column_names)
        index_str = lpad(i, index_length, '0')
        space_padding = total_length - length(index_str) - length(name) - 3 
        formatted_name = "[" * index_str * " " * repeat(" ", space_padding) * name * "]"
        column_format[name] = formatted_name
    end

    return column_format
end

function Pitch_analysis_buffer(time::Float64, Gas_grid::Dict, Dust_grid::Dict)
    column_names = ["theta","sigma_gas","sigma_dust","rho_g","rho_d","vr_g","vphi_g","vr_d","vphi_d"]
    column_dict = create_column_dict(column_names)

    Sigma_g = Gas_grid["Sigma"]
    rho_g = Gas_grid["rho_m"]
    vr_g = Gas_grid["vr"]
    vphi_g = Gas_grid["vphi"]

    Sigma_d = Dust_grid["Sigma"]
    rho_d = Dust_grid["rho_m"]
    vr_d = Dust_grid["vr"]
    vphi_d = Dust_grid["vphi"]

    radius = Sigma_g.axis[1]
    theta = Sigma_g.axis[2]
    return Pitch_analysis_buffer(
    time,
    rho_g,
    rho_d,
    vr_g,
    vr_d,
    vphi_g,
    vphi_d,
    Sigma_g,
    Sigma_d,
    theta,
    radius,
    column_dict
    )
end

function Spiral_analysis_buffer(time::Float64,radius::Float64, Gas_grid::Dict, Dust_grid::Dict)
    column_names = ["theta","sigma_gas","tilt_gas","sigma_dust","tilt_dust"]
    column_dict = create_column_dict(column_names)

    Sigma_g = Gas_grid["Sigma"]
    tilt_g = Gas_grid["tilt"]

    Sigma_d = Dust_grid["Sigma"]
    tilt_d = Dust_grid["tilt"]


    theta = Sigma_g.axis[1]
    return Spiral_analysis_buffer(
    time,
    radius,
    Sigma_g,
    tilt_g,
    Sigma_d,
    tilt_d,
    theta,
    column_dict
    )
end


function self_check_pitch(input :: Pitch_analysis_buffer)
    iaxis = Array{LinRange}(undef,2)
    iaxis[1] = input.radius
    iaxis[2] = input.theta
    for column in keys(input.column_names)
        if typeof(getfield(input,Symbol(column)))<:gridbackend
            axis = getfield(input, Symbol(column)).axis
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
    iaxis = Array{LinRange}(undef,1)
    iaxis[1] = input.theta
    for column in keys(input.column_names)
        if typeof(getfield(input,Symbol(column)))<:gridbackend
            axis = getfield(input, Symbol(column)).axis
            for i in eachindex(iaxis)
                if !(iaxis[i] == axis[i])
                    error("AxisError: Mismatching of axis in $column")
                end
            end
        end
    end
    @info "Axis self-checking of Spiral_analysis success! "
end

function buffer2output(buffer_struct::Any)
    buffer_type = typeof(buffer_struct)
    if (buffer_type <: Pitch_analysis_buffer)
        self_check_pitch(buffer_struct)
    elseif (buffer_type <: Spiral_analysis_buffer)
        self_check_spiral(buffer_struct)
    else
        @warn "The type of struct $buffer_type is not recorded! It may cause some problem."
    end
    buffer_type = typeof(buffer_struct)
    output_type = Symbol(replace(string(buffer_type), "_buffer" => ""))

    fields = fieldnames(buffer_type)
    converted_values = [convert_field(getfield(buffer_struct, f)) for f in fields]

    return eval(output_type)(converted_values...)
end

function Write_pitch_dat(filename::String, data::Pitch_analysis)
    @info "Starting writting the output file for pitch analysis."
    numberfile = extract_number(filename)
    output_filename = "pitch_phjl_" * numberfile * ".dat"
    open(output_filename, "w") do f
        for i in eachindex(data.radius)
            @printf(f,"# Analysis data at t = %15.10E \n",data.time)
            @printf(f,"# theta_rad = %15.10E \n",data.radius[i])
            println(f,"# ",join(values(data.column_names),"   "))
            for j in eachindex(data.theta)
                for key in keys(data.column_names)
                    skey = Symbol(key)
                    if typeof(getfield(data,skey)) <: Vector
                        @printf(f,"%18.10E ",getfield(data,skey)[j])
                    else
                        @printf(f,"%18.10E ",getfield(data,skey)[i,j])
                    end
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
        println(f,"# ",join(values(data.column_names),"   "))
        for j in eachindex(data.theta)
            for key in keys(data.column_names)
                skey = Symbol(key)
                @printf(f,"%18.10E ",getfield(data,skey)[j])
            end
            println(f,"")
        end
    end
    @info "Finished writting the output file for spiral analysis."
end



for name in names(@__MODULE__; all=true)
    if name ∉ (:include, :eval) && isdefined(@__MODULE__, name)
        @eval export $name
    end
end