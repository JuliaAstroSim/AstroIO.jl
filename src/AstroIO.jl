module AstroIO

_precompile_(true)

using Unitful, UnitfulAstro
using FileIO, CSV, HDF5, JLD
using StaticArrays

using PhysicalParticles



include("CSV.jl")
include("Gadget.jl")

end # module
