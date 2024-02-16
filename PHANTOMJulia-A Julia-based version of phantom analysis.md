# PHANTOMJulia-A Julia-based version of phantom analysis

## Introduction

A tool for analysis phantom-based hydrodynamics simulation result with Julia.

For those who does not familiar with Fortran but refuse enduring the lack effencicy of Python.

The file reading system for the raw data is inspired by **sarrecen**

## Structure of the Project

**11 files** in **src/** for every part of analysis. 

1. **module_initialization.jl**: Importing necessary package of Julia for analysis(**DataFrames.jl**,**LinearAlgebra.jl**,...etc.)
2. **phjldataframe.jl**: The raw data storage system(type: *phjlRawDataFrame*) and those functions that can modify the *phjlRawDataFrame* directly.
3. **kernel_function.jl**: Kernel function of SPH.
4. **grid.jl**: The grid generator of analysis(type: *gridbackend*). There has a function *generate_coordinate_grid* that can generate the cooresponding coordinate from *gridbackend* for each point  
5. **mathematical_tools.jl**: Some common mathmatical operation. Such as coordinate transformation, integration...
6. **read_phantom.jl**: File reading system. Basically the same as how **sarrecen** does.
7. **logging.jl**: Some information printing out function.
8. **physical_quantity.jl**: SPH Intepolation for **a single point**. 
9. **Grid-dependent_analysis.jl** and **Ring-dependent_analysis.jl**: Analysis the SPH interpolation for the points on the grid by using **physical_quantity.jl**. The **Grid-dependent** is for the whole disc(*pitch_00XXX.dat*), and the **Ring-dependent** is for a ring with a given radius(*spiral_00XXX.dat*).
10. **Extract_data.jl**: Verifing the data, and writing out the result.

## *struct* in Julia

1. *phjlRawDataFrame*: Similar to class *SarrecenDataFrame* in Python. To make the analysis possible.

   ~~~julia
   struct phjlRawDataFrame
       dfdata :: DataFrame  					  # Table of the result
       params :: Dict				   				# Some of the information of simulation
   end
   ~~~

2. *gridbackend*: A structure that contain the information of the grid. 

   ~~~julia
   struct gridbackend
       grid :: Array         					# The result for each point on the grid (e.g.surface density, density... etc.)
       axis :: Vector{LinRange}				# The axis for each direction of the grid (e.g. radius, angle... etc.)
       dimension :: Vector{Int}				# The dimention of the grid. Should be as same as length(axis)
   end
   
   ~~~

3. *XXX_analysis_buffer*: Used for verifying the axis(All of the axis in output should be the same)

4. *XXX_analysis*: The result struct. Can be saved as txt(phantom analysis format.)



## Usage

There has a example of usage in **Analysis_main.jl**

Step 1. import PHANTOMJulia by 

~~~julia
include("./PHANTOMJulia.jl")
using .PHANTOMJulia
~~~

Step 2. Read file 

~~~julia
pjdf_list = read_phantom(filename, "sinks") #For the seperation of sinks and particles (2 files)
pjdf_list = read_phantom(filename, "all") 	#For the seperation of sinks and different particles (depend on the number of type of the particles)
~~~

Step 3. Change the coordinate(change to primary star based coordinate from COM, rotate the coordinate so that the z-axis would be proportional to the angular momentum of the primary disc), add the necessary quantity for SPH intepolation.

Step 4. Do the analysis(See **Analysis_main.jl**), store in a dict

Step 5. Generate a *XXX_analysis_buffer* structure, verifing the axis.

Step 6. Transfer to the *XXX_analysis* structure and write the output.