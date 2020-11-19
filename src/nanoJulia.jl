module nanoJulia
using Statistics, FASTX, DataFrames, Printf, Formatting, BioAlignments, XAM, Plots, HDF5
export nanoread
# Readinfo types
include("datatype.jl")
include("utilities.jl")
include("record_handler.jl")
include("nanoread.jl")
include("info_handler.jl")
include("statistics_handler.jl")

end
