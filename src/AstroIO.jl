module AstroIO

__precompile__(true)

using Unitful, UnitfulAstro
using FileIO, JLD
using StaticArrays
using IterTools

using PhysicalParticles

import Base: show
import FileIO: load, save

export
    # Base
    show,

    # FileIO
    load, save,

    HeaderGadget2, KeysGadget2,

    read_gadget2, write_gadget2


include("CSV.jl")
include("Gadget.jl")
include("PrettyPrint.jl")

end # module
