"""
    write_csv(filename, particles, units = uAstro)
    write_csv(filename, particles, units = uAstro)   # 3D particles (Star / Star3D-like)
    write_csv(filename, data::Union{Array,StructArray}, units = uAstro)

Write particle data to a CSV file. The output suffix, header columns and
per-row field set depend on the type of `particles`:

- `Star2D` (or any `T <: Star2D`): writes `"$filename.SPHGas2D.csv"`
  with columns `id, x, y, vx, vy, ax, ay, oldacc, m, Ti_endstep,
  Ti_begstep, GravCost, Potential` (lengths/velocities/accelerations are
  converted from `units`).
- 3D `Star` (or any `T <: Star`): writes `"$filename.SPHGas.csv"`
  with the same columns extended to 3 spatial dimensions.
- A generic `Array` / `StructArray` of mixed types: writes `"$filename.csv"`
  using a 5-column (point-only) or 8-column (full particle) header
  depending on the element type.

Pass `units = nothing` to skip unit conversion and store raw numbers.
"""
function write_csv(filename::AbstractString, particles::Union{Array{T,N}, StructArray{T,N,NT,Tu}}, units = uAstro) where T <: Star2D where N where NT where Tu
    uLength = getuLength(units)
    uVel = getuVel(units)
    uAcc = getuAcc(units)
    uMass = getuMass(units)
    uPotential = getuEnergyUnit(units)

    open("$filename.SPHGas2D.csv", "w") do f
        write(f, "#id | x y", axisunit(uLength), " | vx vy", axisunit(uVel), " | ax ay oldacc", axisunit(uAcc),
                 " | m", axisunit(uMass), " | Ti_endstep Ti_begstep GravCost | Potential", axisunit(uPotential), " | \n")
        for p in particles
            buffer = @sprintf(
                    "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%.17g\n",
                    p.ID,
                    ustrip_or_raw(uLength, p.Pos.x),
                    ustrip_or_raw(uLength, p.Pos.y),
                    ustrip_or_raw(uVel, p.Vel.x),
                    ustrip_or_raw(uVel, p.Vel.y),
                    ustrip_or_raw(uAcc, p.Acc.x),
                    ustrip_or_raw(uAcc, p.Acc.y),
                    ustrip_or_raw(uMass, p.Mass),
                    p.Ti_endstep,
                    p.Ti_begstep,
                    p.GravCost,
                    ustrip_or_raw(uPotential, p.Potential),
                )
            write(f, buffer)
        end
    end

    return true
end

function write_csv(filename::AbstractString, particles::Union{Array{T,N}, StructArray{T,N,NT,Tu}}, units = uAstro) where T <: Star where N where NT where Tu
    uLength = getuLength(units)
    uVel = getuVel(units)
    uAcc = getuAcc(units)
    uMass = getuMass(units)
    uPotential = getuEnergyUnit(units)

    open("$filename.SPHGas.csv", "w") do f
        write(f, "#id | x y z", axisunit(uLength), " | vx vy vz", axisunit(uVel), " | ax ay az oldacc", axisunit(uAcc),
                 " | m", axisunit(uMass), " | Ti_endstep Ti_begstep GravCost | Potential", axisunit(uPotential), " | \n")
        for p in particles
            buffer = @sprintf(
                    "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%d,%.17g\n",
                    p.ID,
                    ustrip_or_raw(uLength, p.Pos.x),
                    ustrip_or_raw(uLength, p.Pos.y),
                    ustrip_or_raw(uLength, p.Pos.z),
                    ustrip_or_raw(uVel, p.Vel.x),
                    ustrip_or_raw(uVel, p.Vel.y),
                    ustrip_or_raw(uVel, p.Vel.z),
                    ustrip_or_raw(uAcc, p.Acc.x),
                    ustrip_or_raw(uAcc, p.Acc.y),
                    ustrip_or_raw(uAcc, p.Acc.z),
                    ustrip_or_raw(uMass, p.Mass),
                    p.Ti_endstep,
                    p.Ti_begstep,
                    p.GravCost,
                    ustrip_or_raw(uPotential, p.Potential),
                )
            write(f, buffer)
        end
    end

    return true
end

function write_csv(filename::AbstractString, data::Union{Array,StructArray}, units = uAstro)
    uLength = getuLength(units)
    uVel = getuVel(units)
    uAcc = getuAcc(units)
    uMass = getuMass(units)

    # 依据元素类型选择 3D/2D 列布局；Gadget2Particle 等 3D 粒子不再被误判成 2D
    if typeof(first(data)) <: Union{AbstractPoint3D, AbstractParticle3D}
        open("$filename.csv", "w") do f
            write(f, "#id | x y z", axisunit(uLength), " | vx vy vz", axisunit(uVel), " | ax ay az oldacc", axisunit(uAcc), " | m", axisunit(uMass), "\n")
            for v in values(data)
                for p in v
                    buffer = @sprintf(
                                "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                                p.ID,
                                ustrip_or_raw(uLength, p.Pos.x),
                                ustrip_or_raw(uLength, p.Pos.y),
                                ustrip_or_raw(uLength, p.Pos.z),
                                ustrip_or_raw(uVel, p.Vel.x),
                                ustrip_or_raw(uVel, p.Vel.y),
                                ustrip_or_raw(uVel, p.Vel.z),
                                ustrip_or_raw(uAcc, p.Acc.x),
                                ustrip_or_raw(uAcc, p.Acc.y),
                                ustrip_or_raw(uAcc, p.Acc.z),
                                ustrip_or_raw(uMass, p.Mass),
                            )
                    write(f, buffer)
                end
            end
        end
    else # 2D particles
        open("$filename.csv", "w") do f
            write(f, "#id | x y", axisunit(uLength), " | vx vy", axisunit(uVel), " | ax ay oldacc", axisunit(uAcc), " | m", axisunit(uMass), "\n")
            for v in values(data)
                for p in v
                    buffer = @sprintf(
                                "%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                                p.ID,
                                ustrip_or_raw(uLength, p.Pos.x),
                                ustrip_or_raw(uLength, p.Pos.y),
                                ustrip_or_raw(uVel, p.Vel.x),
                                ustrip_or_raw(uVel, p.Vel.y),
                                ustrip_or_raw(uAcc, p.Acc.x),
                                ustrip_or_raw(uAcc, p.Acc.y),
                                ustrip_or_raw(uMass, p.Mass),
                            )
                    write(f, buffer)
                end
            end
        end
    end

    return true
end
