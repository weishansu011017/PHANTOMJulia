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
            data_dict = read_dict(f, "data_dict",String, Array)
            theta = read(f, "theta")
            radius = read(f, "radius")
            column_names = read_dict(f, "column_names", Int, String)
            data = Pitch_analysis(time, data_dict, theta, radius, column_names)
        elseif struct_type == "Spiral_analysis"
            time = read(f, "time")
            radius = read(f, "radius")
            data_dict = read_dict(f, "data_dict",String, Array)
            theta = read(f, "theta")
            column_names = read_dict(f, "column_names", Int, String)
            data = Spiral_analysis(time, radius, data_dict, theta, column_names)
        else
            error("TypeError: Non-supported type of data.")
        end
    end
    return data
end

function determine_label(z::Array, column_name::String, column_label::String="")
    if column_label != ""
        colarbar_label = column_label
    elseif occursin("sigma",column_name) || occursin("Sigma",column_name)
        z .*= 8.8878E+6
        if occursin("_g",column_name)
            colarbar_label = L"$\Sigma_g$ [g cm$^{-2}$]"
        elseif occursin("_g",column_name)
            colarbar_label = L"$\Sigma_d$ [g cm$^{-2}$]"
        else
            colarbar_label = L"$\Sigma$ [g cm$^{-2}$]"
        end
    elseif occursin("rho",column_name)
        z .*= 5.9407E-7
        if occursin("_g",column_name)
            colarbar_label = L"$\rho_g$ [g cm$^{-3}$]"
        elseif occursin("_g",column_name)
            colarbar_label = L"$\rho_d$ [g cm$^{-3}$]"
        else
            colarbar_label = L"$\rho$ [g cm$^{-3}$]"
        end
    elseif occursin("vr",column_name)
        z .*= 2.978E+6
        if occursin("_g",column_name)
            colarbar_label = L"$\vr_g$ [cm s$^{-1}$]"
        elseif occursin("_g",column_name)
            colarbar_label = L"$\vr_d$ [cm s$^{-1}$]"
        else
            colarbar_label = L"$\vr$ [cm s$^{-1}$]"
        end
    elseif occursin("vphi",column_name)
        z .*= 2.978E+6
        if occursin("_g",column_name)
            colarbar_label = L"$\vphi_g$ [cm s$^{-1}$]"
        elseif occursin("_g",column_name)
            colarbar_label = L"$\vphi_d$ [cm s$^{-1}$]"
        else
            colarbar_label = L"$\vphi$ [cm s$^{-1}$]"
        end
    else
        colarbar_label = column_label
    end
    if occursin("∇",column_name)
        z ./= 1.496E13
        colarbar_label = latexstring(L"$\nabla$",colarbar_label)
    end
    return colarbar_label
end

function polar_plot(data::Pitch_analysis, column_index::Int, Log_mode::Bool=false)
    plot_font = "Computer Modern"
    default(fontfamily=plot_font,linewidth=2, framestyle=:box, label=nothing, grid=false)

    radius_array = data.radius
    theta_array = data.theta
    column_name = data.column_names[column_index]
    z = copy(data.data_dict[column_name])
    colormap_label = determine_label(z,column_name)
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