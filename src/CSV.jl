function write_csv(filename::String, particles::Array{T,N}; mode = "astro") where T<:Star2D where N
    f = open("$filename.Star2D.csv", "w")

    if mode == "astro"
        write(f, "#id | x y [kpc] | vx vy [kpc/Gyr] | ax ay oldacc [kpc/Gyr^2] | m [Msun] | Ti_endstep Ti_begstep GravCost | Potential [Msun*kpc^2/Gyr^2]\n")
        for p in particles
            buffer = @sprintf(
                "%d,%f,%f,%f,%f,%f,%f,%f,%d,%d,%f\n",
                p.ID,
                ustrip(Float64, u"kpc", p.Pos.x),
                ustrip(Float64, u"kpc", p.Pos.y),
                ustrip(Float64, u"kpc/Gyr", p.Vel.x),
                ustrip(Float64, u"kpc/Gyr", p.Vel.y),
                ustrip(Float64, u"kpc/Gyr^2", p.Acc.x),
                ustrip(Float64, u"kpc/Gyr^2", p.Acc.y),
                ustrip(Float64, u"Msun", p.Mass),
                p.Ti_endstep,
                p.Ti_begstep,
                ustrip(Float64, u"Msun*kpc^2/Gyr^2", p.Potential),
            )
            write(f, buffer)
        end
    elseif mode == "si"
        write(f, "#id | x y [m] | vx vy [m/s] | ax ay oldacc [m/s^2] | m [kg] | Ti_endstep Ti_begstep GravCost | Potential [kg*m^2/s^2]")
    elseif mode == "cgs"
        write(f, "#id | x y [cm] | vx vy [cm/s] | ax ay oldacc [cm/s^2] | m [g] | Ti_endstep Ti_begstep GravCost | Potential [g*cm^2/s^2]")
    else
        error("Unsupported unit mode: ", mode, "\n  Try these: astro, si, cgs")
    end
    close(f)
end

function write_csv(filename::String, particles::Array{Star}; mode = "astro")

end

function write_csv(filename::String, particles::Array{SPHGas2D}; mode = "astro")

end

function write_csv(filename::String, particles::Array{SPHGas}; mode = "astro")

end

function write_csv(
        filename::String, data::Dict;
        mode = "astro",
        seperate = false,
    )
    if seperate
        for key in keys(data)
            write_csv("$filename.$key.csv", data[key])
        end
    else
        
    end
end

function read_csv(filename::String)

end