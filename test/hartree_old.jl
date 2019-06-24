function hcycle(n::Integer, u::AbstractFloat, temp::AbstractFloat, ham0::Matrix{T};
                rhostart::Matrix{T}) where T<:Number
    dim = size(ham0, 1)
    diaj = dim+1

    uhf = zeros(T, size(ham0))
    rho = zeros(T, size(ham0))
    ham = copy(ham0)

    etotal = Float64(10e5)
    etotal_old = 0
    dmu = 0

    itr = 1
    while abs(etotal - etotal_old) > ENERGY_TOL
        itr += 1
        eigsys = eigen!(Hermitian(ham, :U))

        # From the occupation calculate new energy
        occ = occupation(n, temp, eigsys.values)
        etotal_old = etotal
        etotal = sum(occ .* eigsys.values)

        # Density matrix (gives microcanonical if temp=0)
        canonicaldensmat!(rho, occ, eigsys.vectors)

        # Linear mixing
        rho .= (1 - DENSITY_MIXING) * rho + DENSITY_MIXING * rhostart
        rhostart = copy(rho)

        # Interaction strength times densities
        rho .*= u

        # We only write on the upper tridiagonal due to Hermiticity
        # Add the Hartree term
        # ↑↑ <↓↓>
        uhf[1:2*diaj:end] .= @view rho[1+diaj:2*diaj:end]

        # ↓↓ <↑↑>
        uhf[1+diaj:2*diaj:end] .= @view rho[1:2*diaj:end]

        ham = ham0 + uhf

        # Constant Hartree term
        # mean( <↑↑> <↓↓> )
        dmu = mean(real.(rho[1:2*diaj:end]) .* real.(rho[1+diaj:2*diaj:end]))/u

        # modulates "global" chem. potential
        ham[1:diaj:end] .-= dmu
    end
    return etotal, dmu, real.(rho[1:diaj:end])./u
end


function dotcycle(n::Integer, u::AbstractFloat, temp::AbstractFloat,
                  ham0::Matrix{T}, dotpos::Integer) where T<:Float64
    dim = size(ham0, 1)
    diaj = dim + 1

    uhf = zeros(T, size(ham0))
    uhf_old = zeros(T, size(ham0))
    ham = copy(ham0)

    etotal = Float64(10e5)
    etotal_old = 0

    up = 2*dotpos - 1
    down = 2*dotpos

    dmu = 0

    itr = 1
    while abs(etotal - etotal_old) > ENERGY_TOL
        eigsys = eigen!(Hermitian(ham, :U))

        # From the occupation calculate new energy
        occ = occupation(n, temp, eigsys.values)
        etotal_old = etotal
        etotal = sum(occ .* eigsys.values)

        # Density matrix (gives microcanonical if temp=0)
        rho = canonicaldensmat(occ, eigsys.vectors)

        # Interaction strength times densities
        rho .*= u

        # We only write on the upper tridiagonal due to Hermiticity
        # Add the Hartree term
        # ↑↑ <↓↓>
        uhf[up, up] = rho[down, down]

        # ↓↓ <↑↑>
        uhf[down, down] = rho[up, up]

        # Linear mixing
        uhf .= (1 - DENSITY_MIXING) * uhf + DENSITY_MIXING * uhf_old
        uhf_old .= uhf

        ham = ham0 + uhf

        # Constant Hartree and Fock term
        # <↑↑> <↓↓>
        dmu = rho[up, up] * rho[down, down]/u

        # modulates "global" chem. potential
        ham[1:diaj:end] .-= dmu/n
        itr += 1
    end
    return etotal, dmu
end
