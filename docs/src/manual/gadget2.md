```@meta
CurrentModule = AstroIO
DocTestSetup = quote
    using AstroIO
    using PhysicalParticles
end
```

# Gadget-2 I/O

The [Gadget-2 user guide](https://wwwmpa.mpa-garching.mpg.de/gadget/users-guide.pdf)
defines two binary formats:

- **Format 1** &mdash; Each block is `<size><data><size>` with a fixed label
  determined by the block order.
- **Format 2** &mdash; Each block is `<skip><NAME><size+8><skip><size><data><size>`
  with explicit 4-character labels (`HEAD`, `POS `, `VEL `, `ID  `, `MASS`, `ACCE`,
  `POT `, `U   `, `RHO `, `HSML`, `ENTR`, `TEMP`).

`AstroIO` reads both formats transparently &mdash; the decision is made by
peeking at the first 4 bytes of the file.

## Units

The default units used by Gadget-2 itself are bundled as `uGadget2`. For
simulations it is far more convenient to convert to `uAstro`:

```@repl gadget2
using PhysicalParticles, AstroIO
uGadget2
uAstro
```

Conversion is automatic &mdash; the second positional argument of `read_gadget2`
is the *target* unit, the optional third is the *file* unit (defaults to
`uGadget2`):

```julia
# default: read uGadget2 from file, convert to uAstro
header, data = read_gadget2("snapshot.gadget2", uAstro);

# explicit file unit (same as above)
header, data = read_gadget2("snapshot.gadget2", uAstro, uGadget2);

# read with no unit conversion
header, data = read_gadget2("snapshot.gadget2", nothing, uGadget2);

# both raw (unitless Numbers)
header, data = read_gadget2("snapshot.gadget2", nothing, nothing);
```

!!! warning "Be careful with unit conversions"
    Mixing `nothing` and a unit for `units` / `fileunits` is a common source
    of dimension errors. Pick one combination and stick with it.

## Writing snapshots

```julia
# Write in default uGadget2 units
write_gadget2("output.Gadget2", header, data)

# Write in uAstro units
write_gadget2("output.Gadget2", header, data, uAstro)
```

### Format selection

`write_gadget2` defaults to format 2; pass `format2 = false` for format 1:

```julia
write_gadget2("output.gadget2", header, data; format2 = false)   # format 1
```

### Acceleration and potential blocks

Gadget-2 snapshots optionally store gravitational acceleration (`ACCE`) and
potential (`POT `) per particle. Enable them with the corresponding keyword
arguments:

```julia
write_gadget2_format2("with_acc_pot.gadget2", header, data; acc = true, pot = true)
```

The reader picks these blocks up automatically whenever they are present in
the file, so no extra configuration is needed on the read side.

## Position-only reads

For visualization or coarse previews the position block is enough &mdash; it
skips header parsing, the velocity/ID/mass blocks, and per-particle unit
conversion. The result is a `StructArray` of `PVector`:

```julia
positions = read_gadget2_pos("snapshot.gadget2", uAstro)
```

If the file is in format 2, only the `POS ` block is scanned; for format 1
the read fast-forwards to the third block.

## Headers

`HeaderGadget2` is a mutable struct that mirrors the binary layout of the
Gadget-2 header (256 bytes). The simplest way to construct one is to pass
a `Vector` of particles &mdash; the per-type counts are computed
automatically:

```julia
header = HeaderGadget2(data)
```

You can also build one explicitly by passing keyword arguments:

```julia
header = HeaderGadget2(;
    Counts      = MVector{6,Int32}([n_gas, n_halo, n_disk, n_bulge, n_star, n_bh]),
    time        = 0.0,
    redshift    = 0.0,
    counts_total = Counts,
    nfiles      = 1,
)
```

`AstroIO` also exposes helpers to interrogate an existing header:

- `read_mass_from_header(h)` &mdash; per-type Boolean / `nothing` mask.
- `read_all_mass_from_header(h)` &mdash; `true` if every non-zero type has its
  mass declared in the header.
- `read_any_mass_from_header(h)` &mdash; `true` if at least one type does.
- `read_gadget2_header(filename)` &mdash; read only the header from a file,
  skipping blocks.

## FileIO integration

`FileIO` provides `load` and `save` that dispatch on the `.gadget2` extension:

```julia
using FileIO

header, data = load("snapshot.gadget2")                     # uAstro, uGadget2
header, data = load("snapshot.gadget2", uAstro, uGadget2)   # explicit

save("copy.gadget2", header, data)
save("copy.gadget2", header, data, uAstro)
```

## Working with JLD2

For checkpoints and cross-language data exchange, JLD2 round-trips the
full `StructArray` and `HeaderGadget2` types:

```julia
write_gadget2_jld("snapshot.jld2", header, data)
header, data = read_gadget2_jld("snapshot.jld2")
```
