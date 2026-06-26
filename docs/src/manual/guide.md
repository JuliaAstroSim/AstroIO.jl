```@meta
CurrentModule = AstroIO
DocTestSetup = quote
    using AstroIO
    using PhysicalParticles
end
```

# Package Guide

A practical tour of the main `AstroIO` entry points. For an exhaustive list of
exports see the [API Reference](../lib/Methods.md).

## Installation

From the Julia REPL, press `]` to enter the Pkg REPL mode and run

```julia
pkg> add AstroIO
```

or add from the git repository

```julia
pkg> add https://github.com/JuliaAstroSim/AstroIO.jl
```

Run the test suite with

```julia
pkg> test AstroIO
```

## Loading the package

```julia
using AstroIO
# optional, needed if you want to construct particles yourself
using PhysicalParticles
```

`AstroIO` re-exports everything from `PhysicalParticles`, so in practice only
`using AstroIO` is required to access both particle types and I/O functions.

## Gadget-2

The file suffixes `gadget2`, `Gadget2`, `GADGET2` are all recognised. Unit
conversion between `uSI`, `uCGS`, `uAstro`, and `uGadget2` is automatic.

```julia
# Read the snapshot, converting file units (uGadget2) to simulation units (uAstro)
header, data = read_gadget2("snapshot.gadget2", uAstro)

# Write a snapshot (preserves the existing header)
write_gadget2("output.Gadget2", header, data)

# If only data is provided, a default header is generated from the particle collections
write_gadget2("output.GADGET2", data)
```

For position-only reads (e.g. for visualization), use the lighter-weight
`read_gadget2_pos`:

```julia
positions = read_gadget2_pos("snapshot.gadget2", uAstro)
```

See [Gadget-2 I/O](gadget2.md) for the full API, including the `format2`,
`acc`, and `pot` keyword arguments.

## FileIO interfaces

[FileIO.jl](https://github.com/JuliaIO/FileIO.jl) provides `load`/`save`
shortcuts that dispatch on the file extension:

```julia
using FileIO

header, data = load("snapshot.gadget2")              # defaults to uAstro / uGadget2
header, data = load("snapshot.gadget2", uAstro, uGadget2)

save("snapshot_copy.gadget2", header, data)
save("snapshot_copy.gadget2", header, data, uAstro)
```

## CSV

`write_csv` produces a `.csv` (or `.SPHGas.csv` / `.SPHGas2D.csv`) file
containing the per-particle data. The file extension is appended automatically,
so the first argument should be a base path without suffix:

```julia
# One star per row, columns: id, x y z, vx vy vz, ax ay az, m, time steps, potential
stars = [Star(uAstro) for _ in 1:10]
write_csv("snapshot", stars, uAstro)              # creates snapshot.SPHGas.csv

# 2D stars -> snapshot.SPHGas2D.csv
stars2d = [Star2D() for _ in 1:10]
write_csv("snapshot2d", stars2d, nothing)
```

A generic overload accepts a heterogeneous `Vector` or `StructArray` of any
`AbstractParticle3D`/`AbstractParticle2D` subtypes and dispatches on the
element type to pick the correct column layout.

## JLD2

JLD2 is convenient for checkpointing because the file format round-trips the
exact `StructArray` and `HeaderGadget2` types:

```julia
# Save a snapshot together with its Gadget-2 header
write_gadget2_jld("snapshot.jld2", header, data)

# Save arbitrary data (no header)
write_jld("plain.jld2", data)

# Reload
header, data = read_gadget2_jld("snapshot.jld2")
data         = read_jld("plain.jld2")
```

## HDF5

The HDF5 backend reads positions of every `PartType*` group in a snapshot:

```julia
pos = read_hdf_pos("snapshot.h5", u"kpc")
```

Writing, header reading, and per-particle data extraction are placeholders
(see [the source](https://github.com/JuliaAstroSim/AstroIO.jl) for the
`HDF5.jl` work-in-progress).

## RAMSES

`write_ramses` dumps a 7-column ASCII file (x y z vx vy vz m) using
`uAstro` lengths / masses and `u"km/s"` velocities:

```julia
write_ramses("ramses.csv", stars)
```

`read_ramses` is currently a stub pending a RAMSES reader.

## Houdini (`.hcsv`)

A single snapshot:

```julia
write_houdini("frame_0001.hcsv", stars, 0.0, uAstro)
```

For appending successive snapshots to the same file use `write_houdini_append`.

To process a series of snapshots at once (e.g. for animation), pass a list of
file indices and a target format:

```julia
using AstroIO: gadget2

write_houdini("movie.hcsv",
              "snapshots/", "snapshot_", collect(0:10:100), ".gadget2",
              gadget2(), uAstro;
              times = collect(0.0:0.01:0.1),
              time_ratio = 1.0,
              pos_ratio = 1.0,
              vel_ratio = 1.0)
```

## ConfParser

`loadconfig` parses an INI file and dispatches on the `format` field to
return either a `IOConfigGadget2` or `IOConfigJLD2`. `loadfromconfig` then
performs the actual read.

```ini
; example.jld2.ini
[io]
filename = ./snapshot.jld2
format   = jld2
label    = data
```

```julia
cfg = loadconfig("example.jld2.ini")
data = loadfromconfig(cfg)
```

## Tools

`renamesuffixs` and `renamereplace` are batch-renaming helpers:

```julia
# Rename all files starting with "test_rename" to have a .ok extension
renamesuffixs("./", "test_rename", ".ok")

# Replace "before" with "after" in every file name containing it
renamereplace("./", "before", "after")
```

See the [Tools page](tools.md) for more details.
