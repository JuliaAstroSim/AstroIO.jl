"""
    read_hdf(filename::AbstractString)

Open `filename` as a read-only HDF5 file and return the underlying
`HDF5.File` handle. The caller is responsible for closing the handle (or
using `do` block syntax). This is the lowest-level entry point —
[`read_hdf_pos`](@ref), [`read_hdf_particles`](@ref) and
[`read_hdf_header`](@ref) build on top of it.
"""
function read_hdf(filename::AbstractString)
    file = h5open(filename, "r")
end

"""
    read_hdf_header(f::IOStream)

!!! warning "Placeholder"
    This function is a stub and currently returns nothing. It is kept as
    part of the public API for forward compatibility with a future header
    reader.
"""
function read_hdf_header(f::IOStream)

end

"""
    read_hdf_particles(f::IOStream)

!!! warning "Placeholder"
    Stub for a future HDF5 particle iterator. Currently returns nothing.
"""
function read_hdf_particles(f::IOStream)

end

"""
    read_hdf_particles(filename::AbstractString, units = uAstro)

!!! warning "Work in progress"
    Only the file is opened; the actual particle iteration is not yet
    implemented. Callers can rely on [`read_hdf_pos`](@ref) instead.
"""
function read_hdf_particles(filename::AbstractString, units = uAstro)
    file = h5open(filename, "r")

    groups = names(file)

    uLength = getuLength(units)

    # Read
end

"""
    read_hdf_pos(filename::AbstractString, u = u"kpc")

Read every `PartType*/Coordinates` dataset in `filename` and return a
flat `Vector{PVector}` of positions. The default length unit `u"kpc"` is
multiplied into each coordinate so the result is unitful.

This is currently the only fully-implemented HDF5 read path; the
remaining [`read_hdf`](@ref) helpers are placeholders.
"""
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

"""
    write_hdf()

!!! warning "Placeholder"
    HDF5 write support is not yet implemented. This stub exists so that
    downstream code can reference the symbol without `MethodError`.
"""
function write_hdf()

end