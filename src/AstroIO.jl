module AstroIO

__precompile__(true)

using Unitful, UnitfulAstro
using FileIO, JLD2
using WriteVTK
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
    read_gadget2_jld, write_gadget2_jld,
    read_jld, write_jld,

    # CSV
    write_csv,
    read_csv

push!(LOAD_PATH, pwd());

include("Gadget.jl")
include("CSV.jl")
include("JLD2.jl")
include("HDF5.jl")
include("PrettyPrint.jl")

end # module
