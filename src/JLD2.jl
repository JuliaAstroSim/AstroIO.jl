"""
    write_gadget2_jld(filename::AbstractString, header::HeaderGadget2, data)

Write a Gadget-2 snapshot (header + particle data) to `filename` using
JLD2. The header and the data are stored under the keys `"header"` and
`"data"` respectively so they can be recovered losslessly with
[`read_gadget2_jld`](@ref).
"""
function write_gadget2_jld(filename::AbstractString, header::HeaderGadget2, data)
    FileIO.save(filename, Dict("header" => header, "data" => data))
    return true
end

"""
    read_gadget2_jld(filename::AbstractString)

Read a JLD2 file previously written by [`write_gadget2_jld`](@ref) and
return the `(header, data)` tuple.
"""
function read_gadget2_jld(filename::AbstractString)
    header, data = FileIO.load(filename, "header", "data")
    return header, data
end

"""
    write_jld(filename::AbstractString, data)

Store `data` in `filename` as a JLD2 file under the key `"data"`. The
companion reader is [`read_jld`](@ref).
"""
function write_jld(filename::AbstractString, data)
    FileIO.save(filename, Dict("data" => data))
    return true
end

"""
    read_jld(filename::AbstractString)

Load a JLD2 file produced by [`write_jld`](@ref) and return the value
previously stored under the key `"data"`.
"""
function read_jld(filename::AbstractString)
    data = FileIO.load(filename, "data")
    return data
end