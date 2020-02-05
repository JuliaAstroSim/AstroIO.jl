function write_csv(filename::String, particles::Array{T,N}, units = uAstro) where T <: Star2D where N
    f = open("$filename.Star2D.csv", "w")

    uLength, uTime, uCurrent, uTemperature, uLuminosity, uMass, uAmount = getunits(units)

    write(f, "#id | x y [$uLength] | vx vy [$uLength/$uTime] | ax ay oldacc [$uLength/$uTime^2] | m [$uMass] | Ti_endstep Ti_begstep GravCost | Potential [$uMass*$uLength^2/$uTime^2]\n")
    for p in particles
        buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                p.ID,
                ustrip(uLength, p.Pos.x),
                ustrip(uLength, p.Pos.y),
                ustrip(uLength / uTime, p.Vel.x),
                ustrip(uLength / uTime, p.Vel.y),
                ustrip(uLength / uTime^2, p.Acc.x),
                ustrip(uLength / uTime^2, p.Acc.y),
                ustrip(uMass, p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(uMass * uLength^2 / uTime^2, p.Potential),
            )
        write(f, buffer)
    end
    
    close(f)
    return true
end

function write_csv(filename::String, particles::Array{T,N}, units = uAstro) where T <: Star where N
    f = open("$filename.Star.csv", "w")

    uLength, uTime, uCurrent, uTemperature, uLuminosity, uMass, uAmount = getunits(units)

    write(f, "#id | x y z [$uLength] | vx vy vz [$uLength/$uTime] | ax ay az oldacc [$uLength/$uTime^2] | m [$uMass] | Ti_endstep Ti_begstep GravCost | Potential [$uMass*$uLength^2/$uTime^2]\n")
    for p in particles
        buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                p.ID,
                ustrip(uLength, p.Pos.x),
                ustrip(uLength, p.Pos.y),
                ustrip(uLength, p.Pos.z),
                ustrip(uLength / uTime, p.Vel.x),
                ustrip(uLength / uTime, p.Vel.y),
                ustrip(uLength / uTime, p.Vel.z),
                ustrip(uLength / uTime^2, p.Acc.x),
                ustrip(uLength / uTime^2, p.Acc.y),
                ustrip(uLength / uTime^2, p.Acc.z),
                ustrip(uMass, p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(uMass * uLength^2 / uTime^2, p.Potential),
            )
        write(f, buffer)
    end
    
    close(f)
    return true
end

function write_csv(filename::String, particles::Array{T,N}, units = uAstro) where T <: SPHGas2D where N
    f = open("$filename.SPHGas2D.csv", "w")

    uLength, uTime, uCurrent, uTemperature, uLuminosity, uMass, uAmount = getunits(units)

    write(f, "#id | x y [$uLength] | vx vy [$uLength/$uTime] | ax ay oldacc [$uLength/$uTime^2] | m [$uMass] | Ti_endstep Ti_begstep GravCost | Potential [$uMass*$uLength^2/$uTime^2] | \n" * 
                 "#Entropy [$uMass*$uLength^2/$uTime^2/$uTemperature] | Density [$uMass/$uLength^2] | Hsml [$uLength] | rvx rvy [$uLength/$uTime] | divv [$uTime^-1] | curlv [$uTime^-1] | dHsmlRho [] | \n" *
                 "#Pressure [$uMass*$uLength^-1*$uTime^-2] | DtEntropy [$uMass*$uLength^3/$uTime^2/$uTemperature] | MaxSignalVel [$uLength/$uTime] |\n")
    for p in particles
        buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                p.ID,
                ustrip(uLength, p.Pos.x),
                ustrip(uLength, p.Pos.y),
                ustrip(uLength / uTime, p.Vel.x),
                ustrip(uLength / uTime, p.Vel.y),
                ustrip(uLength / uTime^2, p.Acc.x),
                ustrip(uLength / uTime^2, p.Acc.y),
                ustrip(uMass, p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(uMass * uLength^2 / uTime^2, p.Potential),

                ustrip(uMass * uLength^2 / uTime^2 / uTemperature, p.Entropy),
                ustrip(uMass / uLength^2, p.Density),
                ustrip(uLength, p.Hsml),
                ustrip(uLength / uTime, p.RotVel.x),
                ustrip(uLength / uTime, p.RotVel.y),
                ustrip(uTime^-1, p.DivVel),
                ustrip(uTime^-1, p.CurlVel),
                p.dHsmlRho,
                ustrip(uMass * uLength^-1 * uTime^-2, p.Pressure),
                ustrip(uMass * uLength^2 / uTime^3 / uTemperature, p.DtEntropy),
                ustrip(uLength / uTime, p.MaxSignalVel)
            )
        write(f, buffer)
    end
    
    close(f)
    return true
end

function write_csv(filename::String, particles::Array{T,N}, units = uAstro) where T <: SPHGas where N
    f = open("$filename.SPHGas.csv", "w")

    uLength, uTime, uCurrent, uTemperature, uLuminosity, uMass, uAmount = getunits(units)

    write(f, "#id | x y z [$uLength] | vx vy vz [$uLength/$uTime] | ax ay az oldacc [$uLength/$uTime^2] | m [$uMass] | Ti_endstep Ti_begstep GravCost | Potential [$uMass*$uLength^2/$uTime^2] | \n" * 
                 "#Entropy [$uMass*$uLength^2/$uTime^2/$uTemperature] | Density [$uMass/$uLength^3] | Hsml [$uLength] | rvx rvy rvz [$uLength/$uTime] | divv [$uTime^-1] | curlv [$uTime^-1] | dHsmlRho [] | \n" *
                 "#Pressure [$uMass*$uLength^-1*$uTime^-2] | DtEntropy [$uMass*$uLength^3/$uTime^2/$uTemperature] | MaxSignalVel [$uLength/$uTime] |\n")
    for p in particles
        buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                p.ID,
                ustrip(uLength, p.Pos.x),
                ustrip(uLength, p.Pos.y),
                ustrip(uLength, p.Pos.z),
                ustrip(uLength / uTime, p.Vel.x),
                ustrip(uLength / uTime, p.Vel.y),
                ustrip(uLength / uTime, p.Vel.z),
                ustrip(uLength / uTime^2, p.Acc.x),
                ustrip(uLength / uTime^2, p.Acc.y),
                ustrip(uLength / uTime^2, p.Acc.z),
                ustrip(uMass, p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(uMass * uLength^2 / uTime^2, p.Potential),

                ustrip(uMass * uLength^2 / uTime^2 / uTemperature, p.Entropy),
                ustrip(uMass / uLength^3, p.Density),
                ustrip(uLength, p.Hsml),
                ustrip(uLength / uTime, p.RotVel.x),
                ustrip(uLength / uTime, p.RotVel.y),
                ustrip(uLength / uTime, p.RotVel.z),
                ustrip(uTime^-1, p.DivVel),
                ustrip(uTime^-1, p.CurlVel),
                p.dHsmlRho,
                ustrip(uMass * uLength^-1 * uTime^-2, p.Pressure),
                ustrip(uMass * uLength^2 / uTime^3 / uTemperature, p.DtEntropy),
                ustrip(uLength / uTime, p.MaxSignalVel)
            )
        write(f, buffer)
    end
    
    close(f)
    return true
end

function write_csv(filename::String, data::Array, units = uAstro)

    uLength, uTime, uCurrent, uTemperature, uLuminosity, uMass, uAmount = getunits(units)

    if typeof(first(data)) <: AbstractPoint3D
        f = open("$filename.csv", "w")

        write(f, "#id | x y z [$uLength] | vx vy vz [$uLength/$uTime] | ax ay az oldacc [$uLength/$uTime^2] | m [$uMass]\n")
        for v in values(data)
            for p in v
                buffer = @sprintf(
                            "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                            p.ID,
                            ustrip(uLength, p.Pos.x),
                            ustrip(uLength, p.Pos.y),
                            ustrip(uLength, p.Pos.z),
                            ustrip(uLength / uTime, p.Vel.x),
                            ustrip(uLength / uTime, p.Vel.y),
                            ustrip(uLength / uTime, p.Vel.z),
                            ustrip(uLength / uTime^2, p.Acc.x),
                            ustrip(uLength / uTime^2, p.Acc.y),
                            ustrip(uLength / uTime^2, p.Acc.z),
                            ustrip(uMass, p.Mass),
                        )
                write(f, buffer)
            end
        end
            
        close(f)
    else # 2D particles
        f = open("$filename.csv", "w")

        write(f, "#id | x y [$uLength] | vx vy [$uLength/$uTime] | ax ay oldacc [$uLength/$uTime^2] | m [$uMass]\n")
        for v in values(data)
            for p in v
                buffer = @sprintf(
                            "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                            p.ID,
                            ustrip(uLength, p.Pos.x),
                            ustrip(uLength, p.Pos.y),
                            ustrip(uLength / uTime, p.Vel.x),
                            ustrip(uLength / uTime, p.Vel.y),
                            ustrip(uLength / uTime^2, p.Acc.x),
                            ustrip(uLength / uTime^2, p.Acc.y),
                            ustrip(uMass, p.Mass),
                            p.Ti_endstep,
                            p.Ti_begstep,
                            ustrip(uMass * uLength^2 / uTime^2, p.Potential),
                        )
                write(f, buffer)
            end
        end
            
        close(f)
    end
        
    @info "Data saved to $filename.csv"
    return true
end