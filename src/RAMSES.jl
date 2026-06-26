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
    uVel = u"km/s"
    uMass = getuMass(units)

    f = open(filename, "w")
    for p in data
        buffer = @sprintf(
            "%f %f %f %f %f %f %f\n",
            ustrip(uLength, p.Pos.x),
            ustrip(uLength, p.Pos.y),
            ustrip(uLength, p.Pos.z),
            ustrip(uVel, p.Vel.x),
            ustrip(uVel, p.Vel.y),
            ustrip(uVel, p.Vel.z),
            ustrip(uMass, p.Mass),
        )
        write(f, buffer)
    end
    close(f)
    return true
end