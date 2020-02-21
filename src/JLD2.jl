function write_gadget2_jld(filename::String, header::HeaderGadget2, data::Array)
    save(filename, Dict("header" => header, "data" => data))
    return true
end

function read_gadget2_jld(filename::String)
    header, data = load(filename, "header", "data")
    return header, data
end

function write_jld(filename::String, data::Array)
    save(filename, Dict("data" => data))
    return true
end

function read_jld(filename::String)
    data = load(filename, "data")
    return data
end

"""
include("AstroIO.jl/src/AstroIO.jl")
using .AstroIO

h,d = read_gadget2("AstroIO.jl/test/gassphere_littleendian.g2")

write_jld("temp/jldtest.jld2", d)


data = load("temp/jldtest.jld2", "data")

@load "temp/jldtest.jld2"


write_jld("temp/jldtest.jld2", d)


data = load("temp/jldtest.jld2", "data")

@load "temp/jldtest.jld2"
"""