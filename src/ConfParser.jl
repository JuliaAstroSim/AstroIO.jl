abstract type IOConfig end

struct IOConfigGadget2 <: IOConfig
    filename::String
    format2::Bool
    iofields::Vector{String}
end

function Base.show(io::IO, ioconfig::IOConfigGadget2)
    print(io,
        """
        Gadget2 I/O config:
            filename: $(ioconfig.filename)
             format2: $(ioconfig.format2)
            iofields: $(ioconfig.iofields)
        """
    )
end

function loadfromconfig(ioconfig::IOConfigGadget2)
    #TODO how to construct Gadget2Block from iofields
    
end



struct IOConfigJLD2 <: IOConfig
    filename::String
    label::String
end

function Base.show(io::IO, ioconfig::IOConfigJLD2)
    print(io,
        """
        JLD2 I/O config:
              filename: $(ioconfig.filename)
            data label: $(ioconfig.label)
        """
    )
end

function loadfromconfig(ioconfig::IOConfigJLD2)
    return FileIO.load(ioconfig.filename, ioconfig.label)
end

function loadconfig(ConfFile::String)
    conf = ConfParse(ConfFile)
    parse_conf!(conf)

    filename = retrieve(conf, "io", "filename")
    format = retrieve(conf, "io", "format")

    if format == "gadget2"
        format2 = retrieve(conf, "io", "format2", Bool)
        iofields = retrieve(conf, "io", "iofields")
        
        # append to 4 Char
        for i in eachindex(iofields)
            Len = length(iofields[i])
            @assert 1 <= Len <= 4 "Field names must have 1~4 characters!"
            iofields[i] *= " "^(4 - Len)
        end

        return IOConfigGadget2(filename, format2, iofields)
    elseif format == "jld2"
        label = retrieve(conf, "io", "label")
        return IOConfigJLD2(filename, label)
    else
        error("Unsupported snapshot format!")
    end
end

