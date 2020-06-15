function write_houdini_header(f::IOStream, pos::Array{T,N}) where T<:AbstractPoint3D where N
    write(f, "id,Px,Py,Pz,time\n")
end

function write_houdini_header(f::IOStream, pos::Array{T,N}) where T<:AbstractParticle3D where N
    write(f, "id,Px,Py,Pz,Vx,Vy,Vz,time\n")
end

function write_houdini_data(f::IOStream, pos::Array{T,N}, time::Float64, units = uAstro;
                            pos_ratio = 1.0,
                            ) where T<:AbstractPoint3D where N
    uLength = getuLength(units)
    for p in pos
        buffer = @sprintf(
            "%d,%f,%f,%f,%f\n",
            p.ID,
            ustrip(uLength, p.x) * pos_ratio,
            ustrip(uLength, p.y) * pos_ratio,
            ustrip(uLength, p.z) * pos_ratio,
            time,
        )
        write(f, buffer)
    end
    flush(f)
end

function write_houdini_data(f::IOStream, particles::Array{T,N}, time::Float64, units = uAstro;
                            pos_ratio = 1.0,
                            vel_ratio = 1.0,
                            ) where T<:AbstractParticle3D where N
    uLength = getuLength(units)
    uVel = getuVel(units)

    # Headers
    for p in particles
        buffer = @sprintf(
            "%d,%f,%f,%f,%f,%f,%f,%f\n",
            p.ID,
            ustrip(uLength, p.Pos.x) * pos_ratio,
            ustrip(uLength, p.Pos.y) * pos_ratio,
            ustrip(uLength, p.Pos.z) * pos_ratio,
            ustrip(uVel, p.Vel.x) * vel_ratio,
            ustrip(uVel, p.Vel.y) * vel_ratio,
            ustrip(uVel, p.Vel.z) * vel_ratio,
            time,
        )
        write(f, buffer)
    end
    flush(f)
end

function write_houdini(filename::AbstractString, data, time::Float64, units = uAstro)
    f = open(filename, "w")
    write_houdini_header(f, data)
    status = write_houdini_data(f, data, time, units)
    close(f)
    return status
end

function write_houdini_append(filename::AbstractString, data, time::Float64, units = uAstro)
    f = open(filename, "a")
    status = write_houdini_data(f, data, time, units)
    close(f)
    return status
end

function write_houdini(OutputFile::AbstractString, folder::AbstractString, filenamebase::AbstractString,
                       Counts::Array{Int64,1}, suffix::AbstractString,
                       FileType::AbstractOutputType, units = uAstro;
                       times = Counts,
                       time_ratio = 1.0,
                       pos_ratio = 1.0,
                       vel_ratio = 1.0,)
    output = open(OutputFile, "w")
    write(output, "id,Px,Py,Pz,Vx,Vy,Vz,time\n")

    progress = Progress(length(Counts), "Loading data and precessing: ")
    for i in eachindex(Counts)
        input = joinpath(folder, string(filenamebase, @sprintf("%04d", Counts[i]), suffix))

        if FileType == gadget2()
            header, data = read_gadget2(input)
        elseif FileType == jld2()
            data = read_jld(input)
        end

        write_houdini_data(output, data, times[i] * time_ratio, units; pos_ratio = pos_ratio, vel_ratio = vel_ratio)

        next!(progress, showvalues = [("iter", i), ("time", times[i]), ("file", input)])
    end
    
    close(output)
    println("Data saved to ", OutputFile)
    return true
end

#write_houdini("HoudiniNiagara.hcsv", "output", "snapshot_", collect(0:10:100), ".gadget2", gadget2(), times = collect(0.0:0.01:0.1))