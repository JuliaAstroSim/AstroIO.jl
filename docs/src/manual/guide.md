# Package Guide

## Installation

From the Julia REPL, type `]` to enter the Pkg REPL mode and run
```julia
pkg> add AstroIO
```
or add from git repository
```julia
pkg> add https://github.com/JuliaAstroSim/AstroIO.jl
```

Test the package by
```julia
pkg> test AstroIO
```

## Basic Usage

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
