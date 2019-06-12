function
temperature(H::Matrix{<:Number}, cycle, args...; rhostart::Matrix{<:Number},
            ne::Int64, u::Float64, Trange::StepRangeLen)

    # Containers
    runs = size(Trange, 1)
    etotal = zeros(Float64, runs)
    hartree = similar(etotal)
    rho = zeros(Float64, runs, size(H, 1))

    for (i, T) in enumerate(Trange)
        println("Run: ", i)
        println("ne: ", ne, " u: ", u, " T: ", T)
        etotal[i], hartree[i], rho[i, :] = cycle(ne, u, T, H, args...;
                                                 rhostart=rhostart)
        println("Particle number: ", sum(rho[i, :]))
    end
    return etotal, hartree, rho
end

function
interaction(H::Matrix{<:Number}, cycle, args...; rhostart::Matrix{<:Number},
            ne::Int64, urange::StepRangeLen, T::Float64)

    # Containers
    runs = size(urange, 1)
    etotal = zeros(Float64, runs)
    hartree = similar(etotal)
    rho = zeros(Float64, runs, size(H, 1))

    for (i, u) in enumerate(urange)
        println("Run: ", i)
        println("ne: ", ne, " u: ", u, " T: ", T)
        etotal[i], hartree[i], rho[i, :] = cycle(ne, u, T, H, args...;
                                                 rhostart=rhostart)
        println("Particle number: ", sum(rho[i, :]))
    end
    return etotal, hartree, rho
end

function
electrons(H::Matrix{<:Number}, cycle, args...;
          nerange::StepRangeLen{Integer, Integer}, u::Float64, T::Float64)

    # Containers
    runs = size(nerange, 1)
    etotal = zeros(Float64, runs)
    hartree = similar(etotal)
    rho = zeros(Float64, runs, size(H, 1))

    for (i, ne) in enumerate(nerange)
        println("Run: ", i)
        println("ne: ", ne, " u: ", u, " T: ", T)
        etotal[i], hartree[i] = cycle(ne, u, T, H, args...)
        println("Particle number: ", sum(rho[i, :]))
    end
    return etotal, hartree, rho
end