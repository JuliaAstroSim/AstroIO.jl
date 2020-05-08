function read_hdf(filename::AbstractString)
    file = h5open(filename, "r")
end

function read_hdf_header(f::IOStream)
    
end

function read_hdf_particles(f::IOStream)
    
end

function read_hdf_particles(filename::AbstractString, units = uAstro)
    file = h5open(filename, "r")

    groups = names(file)

    uLength = getuLength(units)

    # Read
end

function read_hdf_pos(filename::AbstractString, u = u"kpc")
    file = h5open(filename, "r")

    groups = names(file)

    pos = Array{PVector,1}()
    for g in groups
        if startswith(g, "PartType")
            append!(pos, pconvert(Array(file[g]["Coordinates"]) .* u))
        end
    end
    return pos
end

function write_hdf()
    
end