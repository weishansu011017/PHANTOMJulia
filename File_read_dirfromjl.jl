include("jl_module/Module_initialization.jl")
include("jl_module/read_phantom.jl")
include("jl_module/kernel_function.jl")
include("jl_module/physical_quantity.jl")
include("jl_module/Grid-dependent_analysis.jl")


function main()
    filepath = "Test_file/disc_00102"
    prdf_list = read_phantom(filepath,"all")   
    COM2primary(prdf_list,prdf_list[end])
    add_rho(prdf_list[1])
    # println(prdf_list[end].params)
    # println(first(prdf_list[end].dfdata,10))
    # println(first(prdf_list[1].dfdata,10))
    # println(density(prdf_list[1],[42.4955,2.4631,-1.50388],M6_spline,"cylin"))
    println(surface_density(prdf_list[1],[42.4955,2.4631],1.5,M6_spline,100,"cylin"))
    Sigma = surface_density_analysis(prdf_list[1], (10.0,100.0,90),250, 1.3, M6_spline,4)
    heatmap(LinRange(0.0,2*pi,20),LinRange(10.0,100.0,45),Sigma,  projection = :polar, color = :cividis)
    
end


main()




