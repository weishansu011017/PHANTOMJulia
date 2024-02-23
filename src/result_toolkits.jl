function read_dict(f, name, K = String, V = Float64)
    dict = Dict{K, V}()
    g = f[name]
    for key in keys(g)
        val = read(g, key)
        if typeof(key) <: K
            dict[key] = val
        elseif (typeof(key)<: String) && (K <: Int)
            dict[parse(K, key)] = val
        else
            dict[convert(K, key)] = val
        end
    end
    return dict
end

function read_H5DF(filepath::String)
    data = nothing
    h5open(filepath, "r") do f
        struct_type = read(f, "struct_type")
        if struct_type == "Pitch_analysis"
            time = read(f, "time")
            data_dict = read_dict(f, "data_dict",Int, Array)
            theta = read(f, "theta")
            radius = read(f, "radius")
            column_names = read_dict(f, "column_names", Int, String)
            data = Pitch_analysis(time, data_dict, theta, radius, column_names, Dict{String,Bool}())
        elseif struct_type == "Spiral_analysis"
            time = read(f, "time")
            radius = read(f, "radius")
            data_dict = read_dict(f, "data_dict",Int, Array)
            theta = read(f, "theta")
            column_names = read_dict(f, "column_names", Int, String)
            data = Spiral_analysis(time, radius, data_dict, theta, column_names, Dict{String,Bool}())
        else
            error("TypeError: Non-supported type of data.")
        end
    end
    return data
end

function transfer_cgs(data :: Abstract_analysis,umass::Float64=1.9891E+33, udist::Float64=1.496E+13, utime::Float64=5.0227287E+06, year::Bool=true,au::Bool=true)
    function replace_grident_exp(latex_str)
        str = latex_str.s
        regex = r"cm\$\^\{(-?\d+)\}"
        m = match(regex, str)
        if m !== nothing
            exponent = parse(Int, m.captures[1]) - 1
            new_exponent_str = "cm\$^{$exponent}"
            new_str = replace(str, m.match => new_exponent_str)
            return LaTeXString(new_str)
        else
            return latex_str
        end
    end
    try
        _ = data.params["_cgs"]
    catch err
        usigma = umass/(udist^2)
        urho = umass/(udist^3)
        uv = udist/utime

        if year
            data.time *= (utime/31536000)
        else
            data.time *= utime
        end

        if !(au)
            data.radius *= udist
        end

        column_unit = Dict{Int,LaTeXString}()
        for key in keys(data.column_names)
            column_name = data.column_names[key]
            if occursin("sigma",column_name) || occursin("Sigma",column_name)
                data.data_dict[key] *= usigma
                if occursin("_g",column_name)
                    column_unit[key] = L"$\Sigma_g$ [g cm$^{-2}$]"
                elseif occursin("_d",column_name)
                    column_unit[key] = L"$\Sigma_d$ [g cm$^{-2}$]"
                else
                    column_unit[key] = L"$\Sigma$ [g cm$^{-2}$]"
                end
            elseif occursin("rho", column_name)
                data.data_dict[key] *= urho
                if occursin("_g",column_name)
                    column_unit[key] = L"$\rho_g$ [g cm$^{-3}$]"
                elseif occursin("_d",column_name)
                    column_unit[key] = L"$\rho_d$ [g cm$^{-3}$]"
                else
                    column_unit[key] = L"$\rho$ [g cm$^{-3}$]"
                end
            elseif occursin("vr",column_name)
                data.data_dict[key] *= uv
                if occursin("_g",column_name)
                    column_unit[key] = L"$v_{r,g}$ [cm s$^{-1}$]"
                elseif occursin("_d",column_name)
                    column_unit[key] = L"$v_{r,d}$ [cm s$^{-1}$]"
                else
                    column_unit[key] = L"$v_r$ [cm s$^{-1}$]"
                end
            elseif occursin("vphi",column_name)
                data.data_dict[key] *= uv
                if occursin("_g",column_name)
                    column_unit[key] = L"$v_{\phi,g}$ [cm s$^{-1}$]"
                elseif occursin("_d",column_name)
                    column_unit[key] = L"$v_{\phi,d}$ [cm s$^{-1}$]"
                else
                    column_unit[key] = L"$v_{\phi}$ [cm s$^{-1}$]"
                end
            else
                column_unit[key] = L""
            end
            if occursin("∇",column_name)
                data.data_dict[key] /= udist
                column_unit[key] = latexstring(L"$\nabla$",replace_grident_exp(column_unit[key]))
            end
            data.params["column_units"] = deepcopy(column_unit)
        end
    end
    data.params["_cgs"] = true
end

function add_more_label(data::Abstract_analysis, column_index::Int, label::LaTeXString)
    try
        _ = data.params["column_units"][column_index]

        if (data.params["column_units"][column_index] != "")
            while true
                println("Old column unit $(column_unit[column_index]) has found. Are you sure you want to replace it? [y/n]")
                yn = String(readline())
                if yn == "y"
                    break
                elseif yn == "n"
                    return
                else
                    nothing
                end
            end
        end
    catch err
        nothing
    end
    data.params["column_units"][column_index] = label
end

function combine_array(data::Abstract_analysis, column_index::Vector{Int}, Combine_method::Function)
    if length(column_index) == 1
        return data.data_dict[column_index[1]], false
    elseif length(column_index) == 0
        error("CombineError: Missing column index.")
    else
        array_size = size(data.data_dict[column_index[1]])
        array_buffer = zeros(Float64,array_size...,length(column_index))
        result = zeros(Float64,array_size...)
        theta_array = data.theta
        for (i,val) in enumerate(column_index)
            for idx in CartesianIndices(array_size)
                array_buffer[idx,i] = data.data_dict[val][idx]
            end
        end
        for i in 1:array_size[1]
            for j in 1:array_size[2]
                result[i,j] = Combine_method(array_buffer[i,j,:],theta_array[j])
            end
        end
    end
    return result, true
end

function azimuthal_gradient_unit(array::Array, theta::Float64)
    sinθ = sin(theta)
    cosθ = cos(theta)
    norma = norm(array)
    normalized_array = array./norma
    return -sinθ * normalized_array[1] + cosθ * normalized_array[2]
end


function polar_plot(data::Pitch_analysis, column_index::Array{Int}, Combine_method::Function=azimuthal_gradient_unit, label::LaTeXString=L"", Log_mode::Bool=false)
    plot_font = "Computer Modern"
    default(fontfamily=plot_font,linewidth=2, framestyle=:box, label=nothing, grid=false)

    transfer_cgs(data)
    radius_array = data.radius
    theta_array = data.theta
    z,combine_flags = combine_array(data, column_index, Combine_method)

    if (label == L"") && !(combine_flags)
        colormap_label = data.params["column_units"][column_index[1]]
    else
        colormap_label = label
    end
    if Log_mode
        scale = :log
        #colormap_label = latexstring(L"$\log$ ",colormap_label)
    else
        scale = :linear
    end
    colormap_label = latexstring(L"\\ ",colormap_label)
    p = heatmap(theta_array,radius_array, z; 
                projection=:polar,
                color = cgrad(:inferno,scale=scale),
                size = (1200, 700),
                right_margin = 40Plots.mm,
                left_margin = 40Plots.mm,
                colorbar_title=colormap_label,
                colorbar_titlefontsize=22,
                titlefontsize=22,
                guidefontsize=17,
                tickfontsize=15)
    display(p)
end






















for name in names(PHANTOMJulia; all=true)
    if name ∉ (:include, :eval) && isdefined(PHANTOMJulia, name)
        @eval export $name
    end
end