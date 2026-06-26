```@meta
CurrentModule = AstroIO
DocTestSetup = quote
    using AstroIO
    using PhysicalParticles
end
```

# CSV

`AstroIO.write_csv` writes per-particle rows to a CSV file. The first argument
is a *base* path; the actual file name is built by appending one of:

- `.SPHGas.csv`  for `Star` particles
- `.SPHGas2D.csv` for `Star2D` particles
- `.csv` for a heterogeneous `Vector`/`StructArray` of any
  `AbstractParticle3D` / `AbstractParticle2D`

Unit-aware columns are written in the units requested via the third positional
argument (defaults to `uAstro`).

## 3D particles

```julia
stars = [Star(uAstro) for _ in 1:10]
write_csv("snapshot", stars, uAstro)
# writes snapshot.SPHGas.csv
```

Columns:

| # | id | x y z | vx vy vz | ax ay az | m | Ti_endstep Ti_begstep GravCost | Potential |

## 2D particles

```julia
stars2d = [Star2D() for _ in 1:10]
write_csv("snapshot2d", stars2d, nothing)
# writes snapshot2d.SPHGas2D.csv
```

Columns:

| # | id | x y | vx vy | ax ay | m | Ti_endstep Ti_begstep GravCost | Potential |

## Heterogeneous particle lists

Passing a `Vector` (or `StructArray`) of mixed `AbstractParticle3D` types
selects the 3D or 2D column layout based on the element type of the first
entry:

```julia
data = vcat([Star()    for _ in 1:5],
            [Ball()    for _ in 1:5])
write_csv("mixed", data, nothing)
# writes mixed.csv
```

## Unit handling

When `units === nothing`, the function writes raw `Float64` columns. With
a unit, each column is `ustrip`-ed in the requested unit so that downstream
tools can re-parse with their own unit conventions.
