```@meta
CurrentModule = AstroIO
DocTestSetup = quote
    using AstroIO
    using PhysicalParticles
end
```

# Houdini (`.hcsv`)

Houdini's *Geometry Spreadsheet* format is a CSV-like file with one row per
particle. `AstroIO` writes either 5 columns (for `AbstractPoint3D` like
`Ball`) or 8 columns (for full `AbstractParticle3D` like `Star`).

## Single snapshot

```julia
stars = [Star(uAstro) for _ in 1:5]
write_houdini("frame_0001.hcsv", stars, 0.0, uAstro)
```

The CSV header (`id,Px,Py,Pz,Vx,Vy,Vz,time`) is written automatically.

## Appending

For a sequential simulation, use `write_houdini_append` to add new rows
without rewriting the header:

```julia
write_houdini_append("movie.hcsv", stars, 1.0, uAstro)
```

## Series over time

`write_houdini` also has a batch form that loops over a sequence of
snapshot files. It is the most efficient way to build a single `.hcsv`
animation for a series of Gadget-2 snapshots:

```julia
using AstroIO: gadget2

write_houdini(
    "movie.hcsv",
    "snapshots/", "snapshot_", collect(0:10:100), ".gadget2",
    gadget2(), uAstro;
    times        = collect(0.0:0.01:0.1),  # output times, one per input
    time_ratio   = 1.0,                     # additional scaling on the time column
    pos_ratio    = 1.0,                     # scale positions (e.g. unit conversion)
    vel_ratio    = 1.0,                     # scale velocities
)
```

Supported source formats: `gadget2()` and `jld2()` (see
[`AbstractOutputType`](@ref)).

## Time / position / velocity ratios

`pos_ratio` and `vel_ratio` let you bake unit conversions into the file
without modifying the input data. For example, to convert kpc → Mpc:

```julia
write_houdini("movie_mpc.hcsv", stars, 0.0, uAstro; pos_ratio = 1e-3)
```
