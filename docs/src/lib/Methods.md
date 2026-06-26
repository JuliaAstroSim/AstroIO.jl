```@meta
CurrentModule = AstroIO
DocTestSetup = quote
    using AstroIO
end
```

# API Reference

This page lists every public symbol exported by `AstroIO`. Click a name to jump to
its docstring.

## Index

```@index
Pages = ["Methods.md"]
```

## Output types

```@docs
AstroIO.AbstractOutputType
AstroIO.gadget2
AstroIO.hdf5
AstroIO.jld2
```

## Gadget-2

### Data types

```@docs
AstroIO.HeaderGadget2
AstroIO.Gadget2Particle
AstroIO.GadgetKeys
AstroIO.GadgetTypes
```

### Reading

```@docs
AstroIO.read_gadget2
AstroIO.read_gadget2_pos
AstroIO.read_gadget2_header
AstroIO.read_all_mass_from_header
AstroIO.read_any_mass_from_header
AstroIO.read_mass_from_header
AstroIO.count_gadget_types
```

### Writing

```@docs
AstroIO.write_gadget2
AstroIO.write_gadget2_format2
AstroIO.write_gadget2_jld
AstroIO.read_gadget2_jld
AstroIO.generate_gadget2_header
```

### FileIO integration

```@docs
AstroIO.load(::File{format"Gadget2"})
```

`FileIO.save(::File{format"Gadget2"}, header, data[, units])` is also
registered (with `units = uGadget2` as default). It mirrors
`AstroIO.load` but writes a snapshot. See the manual page
`manual/gadget2.md` for usage examples.

## CSV

```@docs
AstroIO.write_csv
```

## RAMSES

```@docs
AstroIO.write_ramses
AstroIO.read_ramses
```

## JLD2

```@docs
AstroIO.read_jld
AstroIO.write_jld
```

## HDF5

!!! note "Work in progress"
    The HDF5 backend is currently limited to reading positions from snapshot
    files following the standard `PartTypeX/Coordinates` layout. Write paths
    are placeholders.

```@docs
AstroIO.read_hdf
AstroIO.read_hdf_pos
AstroIO.read_hdf_header
AstroIO.read_hdf_particles
AstroIO.write_hdf
```

## Houdini

```@docs
AstroIO.write_houdini
AstroIO.write_houdini_append
```

## ConfParser

```@docs
AstroIO.loadconfig
AstroIO.loadfromconfig
AstroIO.IOConfigGadget2
AstroIO.IOConfigJLD2
```

## Tools

```@docs
AstroIO.renamesuffixs
AstroIO.renamereplace
```
