module AstroIO

using PrecompileTools
using Reexport
using Unitful, UnitfulAstro
using FileIO, JLD2, HDF5
using StaticArrays
using IterTools
using Printf
using Distributed
using ProgressMeter
using BangBang
using StructArrays
using Combinatorics: permutations
using ConfParser

@reexport using PhysicalParticles

import Base: show, write
import Unitful: Units
import FileIO: Stream, File

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
    write_gadget2_format2,

    read_gadget2_pos,

    GadgetKeys,
    GadgetTypes,

    # RAMSES
    read_ramses,
    write_ramses,

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

    # ConfParser
    loadconfig,
    loadfromconfig,

    # Tools
    renamereplace,
    renamesuffixs



"""
    GadgetTypes

Canonical ordering of the six Gadget-2 particle collections:
`GAS`, `HALO`, `DISK`, `BULGE`, `STAR`, `BLACKHOLE`. The integer index of
each element matches the block ordering expected by the Gadget-2 binary
header (`npart[1]` = gas, `npart[2]` = halo, …, `npart[6]` = black hole).
"""
const GadgetTypes = [GAS, HALO, DISK, BULGE, STAR, BLACKHOLE]

"""
    GadgetKeys

Default JLD2 / dictionary keys used by the Gadget-2 I/O helpers when an
ordered collection of particle groups is written as a `Dict`. They mirror
[`GadgetTypes`](@ref) one-to-one (`"gases"` ↔ GAS, `"haloes"` ↔ HALO, …).
"""
const GadgetKeys = ["gases", "haloes", "disks", "bulges", "stars", "blackholes"]


"""
    AbstractOutputType

Abstract supertype for the singleton tags used to dispatch on the source
format of batch I/O helpers (for example
[`write_houdini`](@ref)'s series-of-snapshots overload).

Concrete subtypes are [`gadget2`](@ref), [`hdf5`](@ref) and [`jld2`](@ref).
Use `gadget2()` / `hdf5()` / `jld2()` as a lightweight value rather than the
type itself.
"""
abstract type AbstractOutputType end

"""
    gadget2 <: AbstractOutputType

Sentinel value selecting the Gadget-2 binary format in batch I/O helpers
that accept a source-format argument (for example the multi-snapshot
overload of [`write_houdini`](@ref)).
"""
struct gadget2 <: AbstractOutputType end

"""
    hdf5 <: AbstractOutputType

Sentinel value selecting the HDF5 snapshot format in batch I/O helpers.
"""
struct hdf5 <: AbstractOutputType end

"""
    jld2 <: AbstractOutputType

Sentinel value selecting the JLD2 format in batch I/O helpers.
"""
struct jld2 <: AbstractOutputType end

include("Gadget.jl")
include("CSV.jl")
include("RAMSES.jl")
include("JLD2.jl")
include("HDF5.jl")
include("Houdini.jl")
include("PrettyPrint.jl")
include("Tools.jl")
include("ConfParser.jl")

include("precompile.jl")
end # module
