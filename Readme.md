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

```julia
header, data = read_gadget2("snapshot.g2")
header, data = read_gadget2("snapshot.gadget2")

write_gadget2("output.g2")
```

### Use FileIO interfaces

```julia
header, data = load("snapshot.g2")
header, data = load("snapshot.gadget2)
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

## Contribution

Welcome issues and PRs. Need help for other snapshot formats.