# AstroIO.jl

I/O interface for astrophysical simulation codes

[![codecov](https://codecov.io/gh/JuliaAstroSim/AstroIO.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstroSim/AstroIO.jl)
[![][docs-dev-img]][docs-dev-url]

## Installation

```julia
]add AstroIO
```

or

```julia
using Pkg; Pkg.add("AstroIO")
```

or

```julia
using Pkg; Pkg.add("https://github.com/JuliaAstroSim/AstroIO.jl")
```

To test the Package:
```julia
]test AstroIO
```

## Documentation

- [**Dev**][docs-dev-url] &mdash; *documentation of the in-development version.*

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://juliaastrosim.github.io/AstroIO.jl/dev

For beginners, it is highly recommended to read the [documentation of PhysicalParticles.jl](https://juliaastrosim.github.io/PhysicalParticles.jl/dev/).

## Usage

```julia
using AstroIO
```

### Gadget2

Suffixes `gadget2`, `Gadget2`, `GADGET2` are supported

```julia
header, data = read_gadget2("snapshot.gadget2", uAstro)

write_gadget2("output.Gadget2", header, data)

# If only data provided, a default header would be generated
write_gadget2("output.GADGET2", data)
```

Supported units: `uAstro`, `uGadget2`, `uSI`, `uCSG`.

### Use FileIO interfaces

```julia
header, data = load("snapshot.gadget2")
save("FileIO.gadget2", header, data)
```

### Output CSV

```julia
write_csv("output", csv) # No suffix
```

### Save and load with JLD2

```julia
write_gadget2_jld("output.jld2", header, data, uGadget2)
write_jld("NoHeader.jld2", data)

header, data = read_gadget2_jld("output.jld2")
data = read("NoHeader.jld2")

# Or simply use JLD2 interfaces
@load "NoHeader.jld2"
```

## Package ecosystem

- Basic data structure: [PhysicalParticles.jl](https://github.com/JuliaAstroSim/PhysicalParticles.jl)
- File I/O: [AstroIO.jl](https://github.com/JuliaAstroSim/AstroIO.jl)
- Initial Condition: [AstroIC.jl](https://github.com/JuliaAstroSim/AstroIC.jl)
- Parallelism: [ParallelOperations.jl](https://github.com/JuliaAstroSim/ParallelOperations.jl)
- Trees: [PhysicalTrees.jl](https://github.com/JuliaAstroSim/PhysicalTrees.jl)
- Meshes: [PhysicalMeshes.jl](https://github.com/JuliaAstroSim/PhysicalMeshes.jl)
- Plotting: [AstroPlot.jl](https://github.com/JuliaAstroSim/AstroPlot.jl)
- Simulation: [ISLENT](https://github.com/JuliaAstroSim/ISLENT)

## Contribution

Welcome issues and PRs. Need help for other snapshot formats.