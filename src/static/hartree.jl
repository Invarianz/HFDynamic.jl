"""
    ENERGY_TOL

Energy tolerance for the Hartree Fock cycle. Determines break condition.
"""
const ENERGY_TOL = 1e-8


"""
    DENSITY_MIXING

Mixing factor between old and new density matrix. Crucial for convergence.
"""
const DENSITY_MIXING = 0.80


"""
    hfcycle(n, u, temp, ham0; rhostart)

Perform a Hartree fock cycle for given temperature and interaction strength.

# Arguments
- `n::Integer`: Number of electrons
- `u::AbstractFloat`: Interaction strength
- `temp::AbstractFloat`: Temperature of the system
- `ham0::Matrix{T}`: Real or complex Hamiltonian with alternating ↑↓ spin ordering
- `rhostart::Matrix{T}`: Real or complex Density matrix; is used to break the initial 
    field configuration symmetry.
"""
function hfcycle(n::Integer, u::AbstractFloat, temp::AbstractFloat, ham0::Matrix{T};
                 rhostart::Matrix{T}) where T<:Number
    dim = size(ham0, 1)
    diaj = dim+1

    uhf = zeros(T, size(ham0))
    rhoinit = copy(rhostart)
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
        rho .= (1 - DENSITY_MIXING) * rho + DENSITY_MIXING * rhoinit
        rhoinit = copy(rho)

        # Interaction strength times densities
        rho .*= u

        # We only write on the upper tridiagonal due to Hermiticity
        # Add the Hartree term
        # ↑↑ <↓↓>
        uhf[1:2*diaj:end] .= @view rho[1+diaj:2*diaj:end]

        # ↓↓ <↑↑>
        uhf[1+diaj:2*diaj:end] .= @view rho[1:2*diaj:end]

        # Add the Fock term
        # ↑↓ <↓↑>; ↓↑ <↑↓> is implicit due to Hermiticity
        # Only write upper!! half of uhf (above we use Hermitian(ham, :U))
        uhf[2:2*diaj:end] .= @view rho[diaj:2*diaj:end]

        ham = ham0 + uhf

        # Constant Hartree and Fock term
        # mean( <↑↑> <↓↓> + <↑↓> <↓↑> )
        dmu = mean(real.(rho[1:2*diaj:end]) .* real.(rho[1+diaj:2*diaj:end]) .+
                   real.(rho[2:2*diaj:end]) .* real.(rho[diaj:2*diaj:end]))/u

        # modulates "global" chem. potential
        ham[1:diaj:end] .-= dmu
    end
    return etotal, dmu, real.(rho[1:diaj:end])./u
end