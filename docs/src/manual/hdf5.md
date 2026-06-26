```@meta
CurrentModule = AstroIO
DocTestSetup = quote
    using AstroIO
end
```

# HDF5

!!! warning "Work in progress"
    The HDF5 backend currently only supports reading positions. Other entry
    points (`read_hdf`, `read_hdf_header`, `read_hdf_particles`, `write_hdf`)
    are placeholders awaiting implementation.

## Reading positions

`read_hdf_pos` scans every `PartType*` group in the snapshot and concatenates
the `Coordinates` datasets into a single `Vector{PVector}`:

```julia
pos = read_hdf_pos("snapshot.h5", u"kpc")
```

The result can be plotted directly, used to drive a `MeshCartesianStatic`,
or fed back into `read_gadget2_pos`-style helpers for further processing.

## Roadmap

The remaining HDF5 entry points are stubs that should mirror the Gadget-2
API once implemented:

| Function | Gadget-2 equivalent |
| -------- | ------------------ |
| `read_hdf_header`     | `read_gadget2_header` |
| `read_hdf_particles`  | `read_gadget2` (full) |
| `write_hdf`           | `write_gadget2` |
