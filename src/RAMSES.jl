"""
    read_ramses(filename::String)

!!! warning "Placeholder"
    RAMSES reader is not yet implemented. The signature is reserved so
    downstream code can reference the symbol.
"""
function read_ramses(filename::String)
    #TODO
end

"""
    write_ramses(filename::String, data, units = uAstro)

Write particle data to a RAMSES-style ASCII file. Each line stores
`x y z vx vy vz m` in the chosen units (positions/mass from `units`,
velocities in km/s). Use [`read_ramses`](@ref) to recover the data
once the reader is implemented.
"""
function write_ramses(filename::String, data, units = uAstro)
    uLength = getuLength(units)
    uVel = isnothing(units) ? nothing : u"km/s"
    uMass = getuMass(units)

    open(filename, "w") do f
        for p in data
            buffer = @sprintf(
                "%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
                ustrip_or_raw(uLength, p.Pos.x),
                ustrip_or_raw(uLength, p.Pos.y),
                ustrip_or_raw(uLength, p.Pos.z),
                ustrip_or_raw(uVel, p.Vel.x),
                ustrip_or_raw(uVel, p.Vel.y),
                ustrip_or_raw(uVel, p.Vel.z),
                ustrip_or_raw(uMass, p.Mass),
            )
            write(f, buffer)
        end
    end
    return true
end
