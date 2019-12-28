function write_csv(filename::String, particles::Array{T,N}; mode = "astro") where T<:Star2D where N
    f = open("$filename.Star2D.csv", "w")

    if mode == "astro"
        write(f, "#id | x y [kpc] | vx vy [kpc/Gyr] | ax ay oldacc [kpc/Gyr^2] | m [Msun] | Ti_endstep Ti_begstep GravCost | Potential [Msun*kpc^2/Gyr^2]\n")
        for p in particles
            buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                p.ID,
                ustrip(u"kpc", p.Pos.x),
                ustrip(u"kpc", p.Pos.y),
                ustrip(u"kpc/Gyr", p.Vel.x),
                ustrip(u"kpc/Gyr", p.Vel.y),
                ustrip(u"kpc/Gyr^2", p.Acc.x),
                ustrip(u"kpc/Gyr^2", p.Acc.y),
                ustrip(u"Msun", p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(u"Msun*kpc^2/Gyr^2", p.Potential),
            )
            write(f, buffer)
        end
    elseif mode == "si"
        write(f, "#id | x y [m] | vx vy [m/s] | ax ay oldacc [m/s^2] | m [kg] | Ti_endstep Ti_begstep GravCost | Potential [kg*m^2/s^2]\n")
        for p in particles
            buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                p.ID,
                ustrip(u"m", p.Pos.x),
                ustrip(u"m", p.Pos.y),
                ustrip(u"m/s", p.Vel.x),
                ustrip(u"m/s", p.Vel.y),
                ustrip(u"m/s^2", p.Acc.x),
                ustrip(u"m/s^2", p.Acc.y),
                ustrip(u"kg", p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(u"kg*m^2/s^2", p.Potential),
            )
            write(f, buffer)
        end
    elseif mode == "cgs"
        write(f, "#id | x y [cm] | vx vy [cm/s] | ax ay oldacc [cm/s^2] | m [g] | Ti_endstep Ti_begstep GravCost | Potential [g*cm^2/s^2]\n")
        for p in particles
            buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                p.ID,
                ustrip(u"cm", p.Pos.x),
                ustrip(u"cm", p.Pos.y),
                ustrip(u"cm/s", p.Vel.x),
                ustrip(u"cm/s", p.Vel.y),
                ustrip(u"cm/s^2", p.Acc.x),
                ustrip(u"cm/s^2", p.Acc.y),
                ustrip(u"g", p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(u"g*cm^2/s^2", p.Potential),
            )
            write(f, buffer)
        end
    else
        error("Unsupported unit mode: ", mode, "\n  Try these: astro, si, cgs")
    end
    close(f)
end

function write_csv(filename::String, particles::Array{T,N}; mode = "astro") where T<:Star where N
    f = open("$filename.Star.csv", "w")

    if mode == "astro"
        write(f, "#id | x y z [kpc] | vx vy vz [kpc/Gyr] | ax ay az oldacc [kpc/Gyr^2] | m [Msun] | Ti_endstep Ti_begstep GravCost | Potential [Msun*kpc^2/Gyr^2]\n")
        for p in particles
            buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                p.ID,
                ustrip(u"kpc", p.Pos.x),
                ustrip(u"kpc", p.Pos.y),
                ustrip(u"kpc", p.Pos.z),
                ustrip(u"kpc/Gyr", p.Vel.x),
                ustrip(u"kpc/Gyr", p.Vel.y),
                ustrip(u"kpc/Gyr", p.Vel.z),
                ustrip(u"kpc/Gyr^2", p.Acc.x),
                ustrip(u"kpc/Gyr^2", p.Acc.y),
                ustrip(u"kpc/Gyr^2", p.Acc.z),
                ustrip(u"Msun", p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(u"Msun*kpc^2/Gyr^2", p.Potential),
            )
            write(f, buffer)
        end
    elseif mode == "si"
        write(f, "#id | x y z [m] | vx vy vz [m/s] | ax ay az oldacc [m/s^2] | m [kg] | Ti_endstep Ti_begstep GravCost | Potential [kg*m^2/s^2]\n")
        for p in particles
            buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                p.ID,
                ustrip(u"m", p.Pos.x),
                ustrip(u"m", p.Pos.y),
                ustrip(u"m", p.Pos.z),
                ustrip(u"m/s", p.Vel.x),
                ustrip(u"m/s", p.Vel.y),
                ustrip(u"m/s", p.Vel.z),
                ustrip(u"m/s^2", p.Acc.x),
                ustrip(u"m/s^2", p.Acc.y),
                ustrip(u"m/s^2", p.Acc.z),
                ustrip(u"kg", p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(u"kg*m^2/s^2", p.Potential),
            )
            write(f, buffer)
        end
    elseif mode == "cgs"
        write(f, "#id | x y z [cm] | vx vy vz [cm/s] | ax ay az oldacc [cm/s^2] | m [g] | Ti_endstep Ti_begstep GravCost | Potential [g*cm^2/s^2]\n")
        for p in particles
            buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                p.ID,
                ustrip(u"cm", p.Pos.x),
                ustrip(u"cm", p.Pos.y),
                ustrip(u"cm", p.Pos.z),
                ustrip(u"cm/s", p.Vel.x),
                ustrip(u"cm/s", p.Vel.y),
                ustrip(u"cm/s", p.Vel.z),
                ustrip(u"cm/s^2", p.Acc.x),
                ustrip(u"cm/s^2", p.Acc.y),
                ustrip(u"cm/s^2", p.Acc.z),
                ustrip(u"g", p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(u"g*cm^2/s^2", p.Potential),
            )
            write(f, buffer)
        end
    else
        error("Unsupported unit mode: ", mode, "\n  Try these: astro, si, cgs")
    end
    close(f)
end

function write_csv(filename::String, particles::Array{T, N}; mode = "astro") where T<:SPHGas2D where N
    f = open("$filename.SPHGas2D.csv", "w")

    if mode == "astro"
        write(f, "#id | x y [kpc] | vx vy [kpc/Gyr] | ax ay oldacc [kpc/Gyr^2] | m [Msun] | Ti_endstep Ti_begstep GravCost | Potential [Msun*kpc^2/Gyr^2] | \n" * 
                 "#Entropy [Msun*kpc^2/Gyr^2/K] | Density [Msun/kpc^2] | Hsml [kpc] | rvx rvy [kpc/Gyr] | divv [Gyr^-1] | curlv [Gyr^-1] | dHsmlRho [] | \n" *
                 "#Pressure [Msun*kpc^-1*Gyr^-2] | DtEntropy [Msun*kpc^3/Gyr^2/K] | MaxSignalVel [kpc/Gyr] |\n")
        for p in particles
            buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                p.ID,
                ustrip(u"kpc", p.Pos.x),
                ustrip(u"kpc", p.Pos.y),
                ustrip(u"kpc/Gyr", p.Vel.x),
                ustrip(u"kpc/Gyr", p.Vel.y),
                ustrip(u"kpc/Gyr^2", p.Acc.x),
                ustrip(u"kpc/Gyr^2", p.Acc.y),
                ustrip(u"Msun", p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(u"Msun*kpc^2/Gyr^2", p.Potential),

                ustrip(u"Msun*kpc^2/Gyr^2/K", p.Entropy),
                ustrip(u"Msun/kpc^2", p.Density),
                ustrip(u"kpc", p.Hsml),
                ustrip(u"kpc/Gyr", p.RotVel.x),
                ustrip(u"kpc/Gyr", p.RotVel.y),
                ustrip(u"Gyr^-1", p.DivVel),
                ustrip(u"Gyr^-1", p.CurlVel),
                p.dHsmlRho,
                ustrip(u"Msun*kpc^-1*Gyr^-2", p.Pressure),
                ustrip(u"Msun*kpc^2/Gyr^3/K", p.DtEntropy),
                ustrip(u"kpc/Gyr", p.MaxSignalVel)
            )
            write(f, buffer)
        end
    elseif mode == "si"
        write(f, "#id | x y [m] | vx vy [m/s] | ax ay oldacc [m/s^2] | m [kg] | Ti_endstep Ti_begstep GravCost | Potential [kg*m^2/s^2] | \n" *
                 "#Entropy [kg*m^2/s^2/K] | Density [kg/m^2] | Hsml [m] | rvx rvy [m/s] | divv [s^-1] | curlv [s^-1] | dHsmlRho [] | \n" *
                 "#Pressure [kg*m^-1*s^-2] | DtEntropy [kg*m^2/s^3/K] | MaxSignalVel [m/s] |\n")
        for p in particles
            buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                p.ID,
                ustrip(u"m", p.Pos.x),
                ustrip(u"m", p.Pos.y),
                ustrip(u"m/s", p.Vel.x),
                ustrip(u"m/s", p.Vel.y),
                ustrip(u"m/s^2", p.Acc.x),
                ustrip(u"m/s^2", p.Acc.y),
                ustrip(u"kg", p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(u"kg*m^2/s^2", p.Potential),

                ustrip(u"kg*m^2/s^2/K", p.Entropy),
                ustrip(u"kg/m^2", p.Density),
                ustrip(u"m", p.Hsml),
                ustrip(u"m/s", p.RotVel.x),
                ustrip(u"m/s", p.RotVel.y),
                ustrip(u"s^-1", p.DivVel),
                ustrip(u"s^-1", p.CurlVel),
                dHsmlRho,
                ustrip(u"kg*m^-1*s^-2", p.Pressure),
                ustrip(u"kg*m^2/s^3/K", p.DtEntropy),
                ustrip(u"m/s", p.MaxSignalVel)
            )
            write(f, buffer)
        end
    elseif mode == "cgs"
        write(f, "#id | x y [cm] | vx vy [cm/s] | ax ay oldacc [cm/s^2] | m [g] | Ti_endstep Ti_begstep GravCost | Potential [g*cm^2/s^2] | \n" *
                 "#Entropy [g*cm^2/s^2/K] | Density [g/cm^2] | Hsml [cm] | rvx rvy [cm/s] | divv [s^-1] | curlv [s^-1] | dHsmlRho [] | \n" *
                 "#Pressure [g*cm^-1*s^-2] | DtEntropy [g*cm^2/s^3/K] | MaxSignalVel [cm/s] |\n")
        for p in particles
            buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                p.ID,
                ustrip(u"cm", p.Pos.x),
                ustrip(u"cm", p.Pos.y),
                ustrip(u"cm/s", p.Vel.x),
                ustrip(u"cm/s", p.Vel.y),
                ustrip(u"cm/s^2", p.Acc.x),
                ustrip(u"cm/s^2", p.Acc.y),
                ustrip(u"g", p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(u"g*cm^2/s^2", p.Potential),

                ustrip(u"g*cm^2/s^2/K", p.Entropy),
                ustrip(u"g/cm^2", p.Density),
                ustrip(u"cm", p.Hsml),
                ustrip(u"cm/s", p.RotVel.x),
                ustrip(u"cm/s", p.RotVel.y),
                ustrip(u"s^-1", p.DivVel),
                ustrip(u"s^-1", p.CurlVel),
                dHsmlRho,
                ustrip(u"g*cm^-1*s^-2", p.Pressure),
                ustrip(u"g*cm^2/s^3/K", p.DtEntropy),
                ustrip(u"cm/s", p.MaxSignalVel)
            )
            write(f, buffer)
        end
    else
        error("Unsupported unit mode: ", mode, "\n  Try these: astro, si, cgs")
    end
    close(f)
end

function write_csv(filename::String, particles::Array{T, N}; mode = "astro") where T<:SPHGas where N
    f = open("$filename.SPHGas.csv", "w")

    if mode == "astro"
        write(f, "#id | x y z [kpc] | vx vy vz [kpc/Gyr] | ax ay az oldacc [kpc/Gyr^2] | m [Msun] | Ti_endstep Ti_begstep GravCost | Potential [Msun*kpc^2/Gyr^2] | \n" * 
                 "#Entropy [Msun*kpc^2/Gyr^2/K] | Density [Msun/kpc^3] | Hsml [kpc] | rvx rvy rvz [kpc/Gyr] | divv [Gyr^-1] | curlv [Gyr^-1] | dHsmlRho [] | \n" *
                 "#Pressure [Msun*kpc^-1*Gyr^-2] | DtEntropy [Msun*kpc^3/Gyr^2/K] | MaxSignalVel [kpc/Gyr] |\n")
        for p in particles
            buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                p.ID,
                ustrip(u"kpc", p.Pos.x),
                ustrip(u"kpc", p.Pos.y),
                ustrip(u"kpc", p.Pos.z),
                ustrip(u"kpc/Gyr", p.Vel.x),
                ustrip(u"kpc/Gyr", p.Vel.y),
                ustrip(u"kpc/Gyr", p.Vel.z),
                ustrip(u"kpc/Gyr^2", p.Acc.x),
                ustrip(u"kpc/Gyr^2", p.Acc.y),
                ustrip(u"kpc/Gyr^2", p.Acc.z),
                ustrip(u"Msun", p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(u"Msun*kpc^2/Gyr^2", p.Potential),

                ustrip(u"Msun*kpc^2/Gyr^2/K", p.Entropy),
                ustrip(u"Msun/kpc^3", p.Density),
                ustrip(u"kpc", p.Hsml),
                ustrip(u"kpc/Gyr", p.RotVel.x),
                ustrip(u"kpc/Gyr", p.RotVel.y),
                ustrip(u"kpc/Gyr", p.RotVel.z),
                ustrip(u"Gyr^-1", p.DivVel),
                ustrip(u"Gyr^-1", p.CurlVel),
                p.dHsmlRho,
                ustrip(u"Msun*kpc^-1*Gyr^-2", p.Pressure),
                ustrip(u"Msun*kpc^2/Gyr^3/K", p.DtEntropy),
                ustrip(u"kpc/Gyr", p.MaxSignalVel)
            )
            write(f, buffer)
        end
    elseif mode == "si"
        write(f, "#id | x y z [m] | vx vy vz [m/s] | ax ay az oldacc [m/s^2] | m [kg] | Ti_endstep Ti_begstep GravCost | Potential [kg*m^2/s^2] | \n" *
                 "#Entropy [kg*m^2/s^2/K] | Density [kg/m^3] | Hsml [m] | rvx rvy rvz [m/s] | divv [s^-1] | curlv [s^-1] | dHsmlRho [] | \n" *
                 "#Pressure [kg*m^-1*s^-2] | DtEntropy [kg*m^2/s^3/K] | MaxSignalVel [m/s] |\n")
        for p in particles
            buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                p.ID,
                ustrip(u"m", p.Pos.x),
                ustrip(u"m", p.Pos.y),
                ustrip(u"m", p.Pos.z),
                ustrip(u"m/s", p.Vel.x),
                ustrip(u"m/s", p.Vel.y),
                ustrip(u"m/s", p.Vel.z),
                ustrip(u"m/s^2", p.Acc.x),
                ustrip(u"m/s^2", p.Acc.y),
                ustrip(u"m/s^2", p.Acc.z),
                ustrip(u"kg", p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(u"kg*m^2/s^2", p.Potential),

                ustrip(u"kg*m^2/s^2/K", p.Entropy),
                ustrip(u"kg/m^3", p.Density),
                ustrip(u"m", p.Hsml),
                ustrip(u"m/s", p.RotVel.x),
                ustrip(u"m/s", p.RotVel.y),
                ustrip(u"m/s", p.RotVel.z),
                ustrip(u"s^-1", p.DivVel),
                ustrip(u"s^-1", p.CurlVel),
                dHsmlRho,
                ustrip(u"kg*m^-1*s^-2", p.Pressure),
                ustrip(u"kg*m^2/s^3/K", p.DtEntropy),
                ustrip(u"m/s", p.MaxSignalVel)
            )
            write(f, buffer)
        end
    elseif mode == "cgs"
        write(f, "#id | x y z [cm] | vx vy vz [cm/s] | ax ay az oldacc [cm/s^2] | m [g] | Ti_endstep Ti_begstep GravCost | Potential [g*cm^2/s^2] | \n" *
                 "#Entropy [g*cm^2/s^2/K] | Density [g/cm^3] | Hsml [cm] | rvx rvy rvz [cm/s] | divv [s^-1] | curlv [s^-1] | dHsmlRho [] | \n" *
                 "#Pressure [g*cm^-1*s^-2] | DtEntropy [g*cm^2/s^3/K] | MaxSignalVel [cm/s] |\n")
        for p in particles
            buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
                p.ID,
                ustrip(u"cm", p.Pos.x),
                ustrip(u"cm", p.Pos.y),
                ustrip(u"cm", p.Pos.z),
                ustrip(u"cm/s", p.Vel.x),
                ustrip(u"cm/s", p.Vel.y),
                ustrip(u"cm/s", p.Vel.z),
                ustrip(u"cm/s^2", p.Acc.x),
                ustrip(u"cm/s^2", p.Acc.y),
                ustrip(u"cm/s^2", p.Acc.z),
                ustrip(u"g", p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(u"g*cm^2/s^2", p.Potential),

                ustrip(u"g*cm^2/s^2/K", p.Entropy),
                ustrip(u"g/cm^3", p.Density),
                ustrip(u"cm", p.Hsml),
                ustrip(u"cm/s", p.RotVel.x),
                ustrip(u"cm/s", p.RotVel.y),
                ustrip(u"cm/s", p.RotVel.z),
                ustrip(u"s^-1", p.DivVel),
                ustrip(u"s^-1", p.CurlVel),
                dHsmlRho,
                ustrip(u"g*cm^-1*s^-2", p.Pressure),
                ustrip(u"g*cm^2/s^3/K", p.DtEntropy),
                ustrip(u"cm/s", p.MaxSignalVel)
            )
            write(f, buffer)
        end
    else
        error("Unsupported unit mode: ", mode, "\n  Try these: astro, si, cgs")
    end
    close(f)
end

function write_csv(
        filename::String, data::Dict;
        mode = "astro",
        seperate = false,
    )
    if seperate
        for key in keys(data)
            write_csv(filename, data[key])
            @info "$key saved to $filename.$key.csv"
        end
    else
        if typeof(first(data)[2]) <: AbstractPoint3D
            f = open("$filename.csv", "w")

            if mode == "astro"
                write(f, "#id | x y z [kpc] | vx vy vz [kpc/Gyr] | ax ay az oldacc [kpc/Gyr^2] | m [Msun] | Ti_endstep Ti_begstep GravCost | Potential [Msun*kpc^2/Gyr^2]\n")
                for v in values(data)
                    for p in v
                        buffer = @sprintf(
                            "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                            p.ID,
                            ustrip(u"kpc", p.Pos.x),
                            ustrip(u"kpc", p.Pos.y),
                            ustrip(u"kpc", p.Pos.z),
                            ustrip(u"kpc/Gyr", p.Vel.x),
                            ustrip(u"kpc/Gyr", p.Vel.y),
                            ustrip(u"kpc/Gyr", p.Vel.z),
                            ustrip(u"kpc/Gyr^2", p.Acc.x),
                            ustrip(u"kpc/Gyr^2", p.Acc.y),
                            ustrip(u"kpc/Gyr^2", p.Acc.z),
                            ustrip(u"Msun", p.Mass),
                            p.Ti_endstep,
                            p.Ti_begstep,
                            ustrip(u"Msun*kpc^2/Gyr^2", p.Potential),
                        )
                        write(f, buffer)
                    end
                end
            elseif mode == "si"
                write(f, "#id | x y z [m] | vx vy vz [m/s] | ax ay az oldacc [m/s^2] | m [kg] | Ti_endstep Ti_begstep GravCost | Potential [kg*m^2/s^2]\n")
                for v in values(data)
                    for p in v
                        buffer = @sprintf(
                            "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                            p.ID,
                            ustrip(u"m", p.Pos.x),
                            ustrip(u"m", p.Pos.y),
                            ustrip(u"m", p.Pos.z),
                            ustrip(u"m/s", p.Vel.x),
                            ustrip(u"m/s", p.Vel.y),
                            ustrip(u"m/s", p.Vel.z),
                            ustrip(u"m/s^2", p.Acc.x),
                            ustrip(u"m/s^2", p.Acc.y),
                            ustrip(u"m/s^2", p.Acc.z),
                            ustrip(u"kg", p.Mass),
                            p.Ti_endstep,
                            p.Ti_begstep,
                            ustrip(u"kg*m^2/s^2", p.Potential),
                        )
                        write(f, buffer)
                    end
                end
            elseif mode == "cgs"
                write(f, "#id | x y z [cm] | vx vy vz [cm/s] | ax ay az oldacc [cm/s^2] | m [g] | Ti_endstep Ti_begstep GravCost | Potential [g*cm^2/s^2]\n")
                for v in values(data)
                    for p in v
                        buffer = @sprintf(
                            "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                            p.ID,
                            ustrip(u"cm", p.Pos.x),
                            ustrip(u"cm", p.Pos.y),
                            ustrip(u"cm", p.Pos.z),
                            ustrip(u"cm/s", p.Vel.x),
                            ustrip(u"cm/s", p.Vel.y),
                            ustrip(u"cm/s", p.Vel.z),
                            ustrip(u"cm/s^2", p.Acc.x),
                            ustrip(u"cm/s^2", p.Acc.y),
                            ustrip(u"cm/s^2", p.Acc.z),
                            ustrip(u"g", p.Mass),
                            p.Ti_endstep,
                            p.Ti_begstep,
                            ustrip(u"g*cm^2/s^2", p.Potential),
                        )
                        write(f, buffer)
                    end
                end
            else
                error("Unsupported unit mode: ", mode, "\n  Try these: astro, si, cgs")
            end
            close(f)
        else # 2D particles
            f = open("$filename.csv", "w")

            if mode == "astro"
                write(f, "#id | x y [kpc] | vx vy [kpc/Gyr] | ax ay oldacc [kpc/Gyr^2] | m [Msun] | Ti_endstep Ti_begstep GravCost | Potential [Msun*kpc^2/Gyr^2]\n")
                for v in values(data)
                    for p in v
                        buffer = @sprintf(
                            "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                            p.ID,
                            ustrip(u"kpc", p.Pos.x),
                            ustrip(u"kpc", p.Pos.y),
                            ustrip(u"kpc/Gyr", p.Vel.x),
                            ustrip(u"kpc/Gyr", p.Vel.y),
                            ustrip(u"kpc/Gyr^2", p.Acc.x),
                            ustrip(u"kpc/Gyr^2", p.Acc.y),
                            ustrip(u"Msun", p.Mass),
                            p.Ti_endstep,
                            p.Ti_begstep,
                            ustrip(u"Msun*kpc^2/Gyr^2", p.Potential),
                        )
                        write(f, buffer)
                    end
                end
            elseif mode == "si"
                write(f, "#id | x y [m] | vx vy [m/s] | ax ay oldacc [m/s^2] | m [kg] | Ti_endstep Ti_begstep GravCost | Potential [kg*m^2/s^2]\n")
                for v in values(data)
                    for p in v
                        buffer = @sprintf(
                            "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                            p.ID,
                            ustrip(u"m", p.Pos.x),
                            ustrip(u"m", p.Pos.y),
                            ustrip(u"m/s", p.Vel.x),
                            ustrip(u"m/s", p.Vel.y),
                            ustrip(u"m/s^2", p.Acc.x),
                            ustrip(u"m/s^2", p.Acc.y),
                            ustrip(u"kg", p.Mass),
                            p.Ti_endstep,
                            p.Ti_begstep,
                            ustrip(u"kg*m^2/s^2", p.Potential),
                        )
                        write(f, buffer)
                    end
                end
            elseif mode == "cgs"
                write(f, "#id | x y [cm] | vx vy [cm/s] | ax ay oldacc [cm/s^2] | m [g] | Ti_endstep Ti_begstep GravCost | Potential [g*cm^2/s^2]\n")
                for v in values(data)
                    for p in particles
                        buffer = @sprintf(
                            "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                            p.ID,
                            ustrip(u"cm", p.Pos.x),
                            ustrip(u"cm", p.Pos.y),
                            ustrip(u"cm/s", p.Vel.x),
                            ustrip(u"cm/s", p.Vel.y),
                            ustrip(u"cm/s^2", p.Acc.x),
                            ustrip(u"cm/s^2", p.Acc.y),
                            ustrip(u"g", p.Mass),
                            p.Ti_endstep,
                            p.Ti_begstep,
                            ustrip(u"g*cm^2/s^2", p.Potential),
                        )
                        write(f, buffer)
                    end
                end
            else
                error("Unsupported unit mode: ", mode, "\n  Try these: astro, si, cgs")
            end
            close(f)
        end
        
        @info "Data saved to $filename.csv"
    end
    return true
end

function read_csv(filename::String)

end