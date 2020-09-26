module AstroIO

using Unitful, UnitfulAstro
using FileIO, JLD2, HDF5
using StaticArrays
using IterTools
using Printf
using Distributed
using ProgressMeter
using Setfield
using BangBang

using PhysicalParticles

import Base: show, write

export
    # Base
    show,

    AbstractOutputType,
        gadget2,
        hdf5,
        jld2,

    # Gadget2
    HeaderGadget2,
    count_gadget_types,
    generate_gadget2_header,
    read_gadget2, write_gadget2,
    read_gadget2_jld, write_gadget2_jld,

    GadgetKeys,
    GadgetTypes,

    # CSV
    write_csv,

    # JLD2
    read_jld, write_jld,

    # HDF5
    read_hdf, write_hdf,
    read_hdf_header,
    read_hdf_pos,
    read_hdf_particles,

    # Houdini
    write_houdini,

    # Tools
    renamesuffixs



GadgetTypes = [GAS, HALO, DISK, BULGE, STAR, BLACKHOLE]
GadgetKeys = [:gases, :haloes, :disks, :bulges, :stars, :blackholes]


abstract type AbstractOutputType end

struct gadget2 <: AbstractOutputType end
struct hdf5 <: AbstractOutputType end
struct jld2 <: AbstractOutputType end

include("Gadget.jl")
include("CSV.jl")
include("JLD2.jl")
include("HDF5.jl")
include("Houdini.jl")
include("PrettyPrint.jl")
include("Tools.jl")

end # module
