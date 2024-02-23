module PHANTOMJulia
# Include the Julia Module
# With the order of level

#Level 1 (Package)
include("./module_initialization.jl")
#Level 2 (SPH Mathematics)
include("./mathematical_tools.jl")
include("./kernel_function.jl")
#Level 3 (Data Structure)
include("./grid.jl")
include("./phjldataframe.jl")
#Level 4 (Sigal point analysis and File read)
include("./physical_quantity.jl")
include("./read_phantom.jl")
#Level 5 (Analysis)
include("./Grid-dependent_analysis.jl")
include("./Ring-dependent_analysis.jl")
#Level 6 (Extract data and information)
include("./logging.jl")
include("./Extract_data.jl")
include("./result_toolkits.jl")
end