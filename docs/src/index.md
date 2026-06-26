```@meta
CurrentModule = AstroIO
DocTestSetup = quote
    using AstroIO
end
```

# AstroIO.jl

*I/O interfaces for astrophysical simulation codes*

`AstroIO.jl` provides a unified file I/O layer for the [JuliaAstroSim](https://github.com/JuliaAstroSim)
ecosystem. It bridges [PhysicalParticles.jl](https://juliaastrosim.github.io/PhysicalParticles.jl/dev/)
particle types with the most common astrophysical file formats, so simulation data can be read,
written, and round-tripped while preserving unit information and physical attributes.

## Features

- **FileIO.jl integration** &mdash; `load` / `save` dispatch on file extension.
- **Multiple snapshot formats**:
  - [Gadget-2](https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf) (format 1 and format 2, with optional acceleration / potential blocks)
  - JLD2 (with or without Gadget-2 header)
  - CSV (per-particle type: `Star`, `Star2D`, `Ball`, `Ball2D`)
  - HDF5 (positions of every `PartType*` group &mdash; work in progress)
  - RAMSES (CSV-style dump)
- **Visualization exports**:
  - Houdini `.hcsv` (per-snapshot and batch over a series)
- **Configuration-driven I/O** via INI files (`ConfParser.jl` integration).
- **Filesystem helpers** for renaming by suffix / substring.
- **Automatic unit conversion** between `uSI`, `uCGS`, `uAstro`, `uGadget2`, and `nothing`.

## Installation

```julia
julia> ]add AstroIO
```

or via Pkg

```julia
julia> using Pkg; Pkg.add("AstroIO")
```

For the latest development version:

```julia
julia> using Pkg; Pkg.add("https://github.com/JuliaAstroSim/AstroIO.jl")
```

## Quick start

```@repl quickstart
using PhysicalParticles, AstroIO

# Read a Gadget-2 snapshot, converting file units (uGadget2) to simulation units (uAstro)
header, data = read_gadget2(joinpath(pkgdir(AstroIO), "test", "gassphere_littleendian.gadget2"), uAstro)

# Inspect
length(data)        # number of particles
header.npart        # particle counts per type
data.Mass[1]        # first particle mass
data.Pos[1]         # first particle position
```

## Manual Outline

```@contents
Pages = ["manual/guide.md", "manual/gadget2.md", "manual/csv.md",
         "manual/jld2.md", "manual/hdf5.md", "manual/houdini.md",
         "manual/confparser.md", "manual/tools.md"]
Depth = 2
```

## Library

```@contents
Pages = ["lib/Methods.md"]
```

## Related packages

`AstroIO` is part of the JuliaAstroSim ecosystem:

- [PhysicalParticles.jl](https://github.com/JuliaAstroSim/PhysicalParticles.jl) &mdash; particle data structures
- [AstroIC.jl](https://github.com/JuliaAstroSim/AstroIC.jl) &mdash; initial conditions
- [PhysicalMeshes.jl](https://github.com/JuliaAstroSim/PhysicalMeshes.jl) &mdash; mesh data structures
- [PhysicalTrees.jl](https://github.com/JuliaAstroSim/PhysicalTrees.jl) &mdash; tree algorithms
- [ParallelOperations.jl](https://github.com/JuliaAstroSim/ParallelOperations.jl) &mdash; parallelism
- [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl) &mdash; visualization
- [ISLENT](https://github.com/JuliaAstroSim/ISLENT) &mdash; full simulation framework
