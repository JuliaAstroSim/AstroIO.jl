# Gadget2 I/O

The default units in Gadget-2 snapshots are collected as `uGadget2`.
However, it is much more convenient to use `uAstro` in simulations:
```@repl gadget2
using PhysicalParticles, AstroIO
uGadget2
uAstro
```

### Basic usage

Units are automatedly converted during I/O:
```julia
# read default uGadget2 units from file and convert to uAstro
header, data = read_gadget2("snapshot.gadget2", uAstro);

# manually set the units in snapshot file
header, data = read_gadget2("snapshot.gadget2", uAstro, uGadget2);

# read default uGadget2 units from file and convert to unitless
header, data = read_gadget2("snapshot.gadget2", nothing, uGadget2);

# write in default uGadget2 units
write_gadget2("output.Gadget2", header, data)

# write in uAstro units
write_gadget2("output.Gadget2", header, data, uAstro)
# and read in uAstro units
header, data = read_gadget2("snapshot.gadget2", uAstro, uAstro);
```
Be careful with these unit conversions!

In some cases (for example, visualization), it is more efficient to only read position data from file:
```julia
read_gadget2_pos("snapshot.gadget2", uAstro)
read_gadget2_pos("snapshot.gadget2", uAstro, uGadget2)
```

### FileIO interfaces

FileIO provides `load` and `save` interfaces to I/O snapshots:
```julia
header, data = load("snapshot.gadget2")
save("FileIO.gadget2", header, data)

header, data = load("snapshot.gadget2", uAstro, uGadget2)
save("FileIO.gadget2", header, data, uGadget2)
```