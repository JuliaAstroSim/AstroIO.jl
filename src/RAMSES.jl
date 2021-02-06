function read_ramses(filename::String)
    #TODO
end

function write_ramses(filename::String, data, units = uAstro)
    uLength = getuLength(units)
    uVel = u"km/s"
    uMass = getuMass(units)

    f = open(filename, "w")
    for p in Iterators.flatten(values(data))
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