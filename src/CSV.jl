function write_csv(filename::AbstractString, particles::Array{T,N}, units = uAstro) where T <: Star2D where N
    f = open("$filename.Star2D.csv", "w")

    uLength = getuLength(units)
    uVel = getuVel(units)
    uAcc = getuAcc(units)
    uMass = getuMass(units)
    uPotential = getuEnergy(units)

    write(f, "#id | x y", axisunit(uLength), " | vx vy", axisunit(uVel), " | ax ay oldacc", axisunit(uAcc),
             " | m", axisunit(uMass), " | Ti_endstep Ti_begstep GravCost | Potential", axisunit(uPotential), "\n")
    for p in particles
        buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                p.ID,
                ustrip(uLength, p.Pos.x),
                ustrip(uLength, p.Pos.y),
                ustrip(uVel, p.Vel.x),
                ustrip(uVel, p.Vel.y),
                ustrip(uAcc, p.Acc.x),
                ustrip(uAcc, p.Acc.y),
                ustrip(uMass, p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(uPotential, p.Potential),
            )
        write(f, buffer)
    end
    
    close(f)
    return true
end

function write_csv(filename::AbstractString, particles::Array{T,N}, units = uAstro) where T <: Star where N
    f = open("$filename.Star.csv", "w")

    uLength = getuLength(units)
    uVel = getuVel(units)
    uAcc = getuAcc(units)
    uMass = getuMass(units)
    uPotential = getuEnergy(units)

    write(f, "#id | x y z", axisunit(uLength), " | vx vy vz", axisunit(uVel), " | ax ay az oldacc", axisunit(uAcc),
             " | m", axisunit(uMass), " | Ti_endstep Ti_begstep GravCost | Potential", axisunit(uPotential), "\n")
    for p in particles
        buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                p.ID,
                ustrip(uLength, p.Pos.x),
                ustrip(uLength, p.Pos.y),
                ustrip(uLength, p.Pos.z),
                ustrip(uVel, p.Vel.x),
                ustrip(uVel, p.Vel.y),
                ustrip(uVel, p.Vel.z),
                ustrip(uAcc, p.Acc.x),
                ustrip(uAcc, p.Acc.y),
                ustrip(uAcc, p.Acc.z),
                ustrip(uMass, p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(uPotential, p.Potential),
            )
        write(f, buffer)
    end
    
    close(f)
    return true
end

function write_csv(filename::AbstractString, particles::Array{T,N}, units = uAstro) where T <: SPHGas2D where N
    f = open("$filename.SPHGas2D.csv", "w")

    uLength = getuLength(units)
    uVel = getuVel(units)
    uAcc = getuAcc(units)
    uMass = getuMass(units)
    uTime = getuTime(units)
    uPotential = getuEnergy(units)
    uEntropy = getuEntropy(units)
    uDensity2D = getuDensity2D(units)
    uPressure = getuPressure(units)

    if isnothing(units)
        uTimeInv = nothing
        uDtEntropy = nothing
    else
        uTimeInv = uTime^-1
        uDtEntropy = uEntropy / uTime
    end

    write(f, "#id | x y", axisunit(uLength), " | vx vy", axisunit(uVel), " | ax ay oldacc", axisunit(uAcc),
             " | m", axisunit(uMass), " | Ti_endstep Ti_begstep GravCost | Potential", axisunit(uPotential), " | \n",
             "#Entropy", axisunit(uEntropy), " | Density", axisunit(uDensity2D), " | Hsml", axisunit(uLength),
             " | rvx rvy", axisunit(uVel), " | divv", axisunit(uTimeInv), " | curlv", axisunit(uTimeInv), " | dHsmlRho", axisunit(uLength), " | \n",
             "#Pressure", axisunit(uPressure), " | DtEntropy", axisunit(uDtEntropy), " | MaxSignalVel", axisunit(uVel), " |\n")
    for p in particles
        buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                p.ID,
                ustrip(uLength, p.Pos.x),
                ustrip(uLength, p.Pos.y),
                ustrip(uVel, p.Vel.x),
                ustrip(uVel, p.Vel.y),
                ustrip(uAcc, p.Acc.x),
                ustrip(uAcc, p.Acc.y),
                ustrip(uMass, p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(uPotential, p.Potential),

                ustrip(uEntropy, p.Entropy),
                ustrip(uDensity2D, p.Density),
                ustrip(uLength, p.Hsml),
                ustrip(uVel, p.RotVel.x),
                ustrip(uVel, p.RotVel.y),
                ustrip(uTimeInv, p.DivVel),
                ustrip(uTimeInv, p.CurlVel),
                ustrip(uLength, p.dHsmlRho),
                ustrip(uPressure, p.Pressure),
                ustrip(uDtEntropy, p.DtEntropy),
                ustrip(uVel, p.MaxSignalVel)
            )
        write(f, buffer)
    end
    
    close(f)
    return true
end

function write_csv(filename::AbstractString, particles::Array{T,N}, units = uAstro) where T <: SPHGas where N
    f = open("$filename.SPHGas.csv", "w")

    uLength = getuLength(units)
    uVel = getuVel(units)
    uAcc = getuAcc(units)
    uMass = getuMass(units)
    uTime = getuTime(units)
    uPotential = getuEnergy(units)
    uEntropy = getuEntropy(units)
    uDensity = getuDensity(units)
    uPressure = getuPressure(units)

    if isnothing(units)
        uTimeInv = nothing
        uDtEntropy = nothing
    else
        uTimeInv = uTime^-1
        uDtEntropy = uEntropy / uTime
    end

    write(f, "#id | x y z", axisunit(uLength), " | vx vy vz", axisunit(uVel), " | ax ay az oldacc", axisunit(uAcc),
             " | m", axisunit(uMass), " | Ti_endstep Ti_begstep GravCost | Potential", axisunit(uPotential), " | \n",
             "#Entropy", axisunit(uEntropy), " | Density", axisunit(uDensity), " | Hsml", axisunit(uLength),
             " | rvx rvy rvz", axisunit(uVel), " | divv", axisunit(uTimeInv), " | curlv", axisunit(uTimeInv), " | dHsmlRho", axisunit(uLength), " | \n",
             "#Pressure", axisunit(uPressure), " | DtEntropy", axisunit(uDtEntropy), " | MaxSignalVel", axisunit(uVel), " |\n")
    for p in particles
        buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                p.ID,
                ustrip(uLength, p.Pos.x),
                ustrip(uLength, p.Pos.y),
                ustrip(uLength, p.Pos.z),
                ustrip(uVel, p.Vel.x),
                ustrip(uVel, p.Vel.y),
                ustrip(uVel, p.Vel.z),
                ustrip(uAcc, p.Acc.x),
                ustrip(uAcc, p.Acc.y),
                ustrip(uAcc, p.Acc.z),
                ustrip(uMass, p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(uPotential, p.Potential),

                ustrip(uEntropy, p.Entropy),
                ustrip(uDensity, p.Density),
                ustrip(uLength, p.Hsml),
                ustrip(uVel, p.RotVel.x),
                ustrip(uVel, p.RotVel.y),
                ustrip(uVel, p.RotVel.z),
                ustrip(uTimeInv, p.DivVel),
                ustrip(uTimeInv, p.CurlVel),
                ustrip(uLength, p.dHsmlRho),
                ustrip(uPressure, p.Pressure),
                ustrip(uDtEntropy, p.DtEntropy),
                ustrip(uVel, p.MaxSignalVel)
            )
        write(f, buffer)
    end
    
    close(f)
    return true
end

function write_csv(filename::AbstractString, data::Array, units = uAstro)

    uLength = getuLength(units)
    uVel = getuVel(units)
    uAcc = getuAcc(units)
    uMass = getuMass(units)

    if typeof(first(data)) <: AbstractPoint3D
        f = open("$filename.csv", "w")

        write(f, "#id | x y z", axisunit(uLength), " | vx vy vz", axisunit(uVel), " | ax ay az oldacc", axisunit(uAcc), " | m", axisunit(uMass), "\n")
        for v in values(data)
            for p in v
                buffer = @sprintf(
                            "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                            p.ID,
                            ustrip(uLength, p.Pos.x),
                            ustrip(uLength, p.Pos.y),
                            ustrip(uLength, p.Pos.z),
                            ustrip(uVel, p.Vel.x),
                            ustrip(uVel, p.Vel.y),
                            ustrip(uVel, p.Vel.z),
                            ustrip(uAcc, p.Acc.x),
                            ustrip(uAcc, p.Acc.y),
                            ustrip(uAcc, p.Acc.z),
                            ustrip(uMass, p.Mass),
                        )
                write(f, buffer)
            end
        end
            
        close(f)
    else # 2D particles
        f = open("$filename.csv", "w")

        write(f, "#id | x y", axisunit(uLength), " | vx vy", axisunit(uVel), " | ax ay oldacc", axisunit(uAcc), " | m", axisunit(uMass), "\n")
        for v in values(data)
            for p in v
                buffer = @sprintf(
                            "%d,%f,%f,%f,%f,%f,%f,%f\n",
                            p.ID,
                            ustrip(uLength, p.Pos.x),
                            ustrip(uLength, p.Pos.y),
                            ustrip(uVel, p.Vel.x),
                            ustrip(uVel, p.Vel.y),
                            ustrip(uAcc, p.Acc.x),
                            ustrip(uAcc, p.Acc.y),
                            ustrip(uMass, p.Mass),
                        )
                write(f, buffer)
            end
        end
            
        close(f)
    end

    return true
end

function write_csv_seperate(filename::AbstractString, data::Dict, units = uAstro)
    for v in values(data)
        if !isempty(v)
            write_csv(filename, v, units)
        end
    end
end

function write_csv(filename::AbstractString, data::Dict, units = uAstro; seperate::Bool = true)
    if seperate
        a = collect(Iterators.flatten(values(data)))
        write_csv(filename, a, units)
    else
        write_csv_seperate(filename, data, units)
    end
end