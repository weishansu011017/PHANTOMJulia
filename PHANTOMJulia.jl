module PHANTOMJulia
# Include the Julia Module
# With the order of level

#Level 1 (Package)
include("./src/module_initialization.jl")
#Level 2 (SPH Mathematics)
include("./src/mathematical_tools.jl")
include("./src/kernel_function.jl")
#Level 3 (Data Structure)
include("./src/grid.jl")
include("./src/phjldataframe.jl")
#Level 4 (Sigal point analysis and File read)
include("./src/physical_quantity.jl")
include("./src/read_phantom.jl")
#Level 5 (Analysis)
include("./src/Grid-dependent_analysis.jl")
include("./src/Ring-dependent_analysis.jl")
#Level 6 (Extract data and information)
include("./src/logging.jl")
include("./src/Extract_data.jl")
end