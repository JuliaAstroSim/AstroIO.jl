function write_vtk_3d(filename::AbstractString, header::HeaderGadget2, data::Dict, units = uAstro)
    uLength, uTime, uCurrent, uTemperature, uLuminosity, uMass, uAmount = getunits(units)

    # Preparing data
    x = zeros(0)
    y = zeros(0)
    z = zeros(0)

    m = zeros(0)

    for p in Iterators.flatten(values(data))
        append!(x, ustrip(uLength, p.Pos.x))
        append!(y, ustrip(uLength, p.Pos.y))
        append!(z, ustrip(uLength, p.Pos.z))

        append!(m, ustrip(uMass, p.Mass))
    end

    # Output
    vtkfile = vtk_grid(filename, x, y, z)
    @show length(x), length(y), length(z), length(m)
    #vtkfile["Mass"] = m
    return vtk_save(vtkfile)
end

function write_vtk_2d(filename::AbstractString, header::HeaderGadget2, data::Dict, units = uAstro)
    
end

function write_vtk(filename::AbstractString, header::HeaderGadget2, data::Dict, units = uAstro)
    if typeof(first(Iterators.flatten(values(data))).Pos) <: AbstractPoint3D
        write_vtk_3d(filename, header, data, units)
    else
        write_vtk_2d(filename, header, data, units)
    end
end

function read_vtk(filename::AbstractString)

end