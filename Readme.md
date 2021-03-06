# AstroIO.jl

I/O interface for astrophysical simulation codes

[![codecov](https://codecov.io/gh/JuliaAstroSim/AstroIO.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaAstroSim/AstroIO.jl)

## Installation

```julia
]add AstroIO
```
or
```julia
]add https://github.com/JuliaAstroSim/AstroIO.jl
```

## Usage

```julia
using AstroIO
```

### Gadget2

Suffixes `gadget2`, `Gadget2`, `GADGET2` are supported

```julia
header, data = read_gadget2("snapshot.gadget2")

write_gadget2("output.Gadget2", header, data)

# If only data provided, a default header would be generated
write_gadget2("output.GADGET2", data)
```

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
write_gadget2_jld("output.jld2", header, data)
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