module AstroIO

__precompile__(true)

using Unitful, UnitfulAstro
using FileIO, JLD
using StaticArrays
using IterTools
using Printf

using PhysicalParticles

import Base: show
import FileIO: load, save

export
    # Base
    show,

    # FileIO
    load, save,

    # Gadget2
    HeaderGadget2, KeysGadget2,
    read_gadget2, write_gadget2,

    # CSV
    write_csv,
    read_csv


include("Gadget.jl")
include("CSV.jl")
include("PrettyPrint.jl")

end # module
