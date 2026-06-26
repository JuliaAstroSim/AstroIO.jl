```@meta
CurrentModule = AstroIO
DocTestSetup = quote
    using AstroIO
    using PhysicalParticles
end
```

# JLD2

JLD2 round-trips the exact `StructArray` of particles and the `HeaderGadget2`
type, which is convenient for checkpoints and cross-language data exchange.

## Save / load

```julia
using AstroIO

header, data = read_gadget2("snapshot.gadget2", uAstro)

# Save header + data
write_gadget2_jld("snapshot.jld2", header, data)
h2, d2 = read_gadget2_jld("snapshot.jld2")

# Save arbitrary data without a header
write_jld("plain.jld2", data)
d3 = read_jld("plain.jld2")
```

## Underlying keys

`write_gadget2_jld` stores the data under the keys `"header"` and `"data"`,
matching the `read_gadget2_jld` reader. `write_jld` / `read_jld` use the
single key `"data"`. The file can also be opened with `JLD2.@load` directly
if you prefer the macro form:

```julia
using JLD2
@load "snapshot.jld2" header data
```

## When to use JLD2 over Gadget-2 binary

- **Pros** &mdash; lossless round-trip of arbitrary Julia types; small files
  thanks to compression; no manual unit handling required.
- **Cons** &mdash; not a community-standard astrophysics format; not portable
  outside Julia.
