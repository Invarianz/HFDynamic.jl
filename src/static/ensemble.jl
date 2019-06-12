"""
    OCCUPATION_TOL

Tolerance for the number of electrons.
"""
const OCCUPATION_TOL = 1e-8


"""
    occupation(npa, temp, en)

Return the occupation of `npa` number of electrons at temperature `temp`
occupying the energy levels `en`, with the Fermi-Dirac distribution.
The energy levels are assumed to be sorted from smallest to biggest.
"""
function
occupation(npa::Integer, temp::T,
           en::Array{T, 1}) where T<:AbstractFloat

    # Zero temperature solution
    if temp ≤ eps(temp)
        occ = zeros(size(en, 1))
        occ[1:npa] .= 1.0
        return occ
    end

    # Estimate μ as Fermi energy at T=0
    muest = en[npa]

    # Minmum and maximum energies
    bandwidth = abs(en[end] - en[1])
    mumin = en[npa] - bandwidth
    mumax = en[npa] + bandwidth

    #Calculate the occupation
    f(e, mu) = 1/(exp((e-mu)/temp)+1)
    occ = map(e -> f(e, muest), en)
    npathermal = sum(occ)

    while abs(npathermal - npa) ≥ OCCUPATION_TOL
        if npathermal - npa > 0
            mumax = muest
            muest = (muest + mumin)/2
        else
            mumin = muest
            muest = (muest + mumax)/2
        end

        occ = map(e -> f(e, muest), en)
        npathermal = sum(occ)
    end
    return occ
end


"""
    canonicaldensmat(occ, ev)

Return the density matrix of a given eigensystem with a given occupation.
See also: ['occupation'](@ref)
"""
function
canonicaldensmat!(rho::Matrix{T}, occ::Vector{<:AbstractFloat}, ev::Matrix{T}) where T<:Number
    rho .= ev * (adjoint(ev) .* occ)
end