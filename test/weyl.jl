module Wrapper
using LinearAlgebra
using PyPlot

include("tools.jl")
include("hartree.jl")
include("eval.jl")

using .Tools
using .Hartree
using .Eval

function
weyl(X::Integer, Y::Integer, Z::Integer)
    su = ComplexF64[1 0; 0 1]
    sx = ComplexF64[0 1; 1 0]
    sy = ComplexF64[0 -im; +im 0]
    sz = ComplexF64[1 0; 0 -1]

    t = -1.0

    on = 2t .* sz
    hopx = (t/2)im .* sx - (t/2) .* sz
    hopy = (t/2)im .* sy - (t/2) .* sz
    hopz = (t/2) .* sz

    chainx = kron(matdiag(true, X, X), on)
    X > 1 && (chainx .+= kron(matdiag(true, X, X, sc=2), hopx))

    chainhopy = kron(matdiag(true, X, X), hopy)

    slicexy = kron(matdiag(true, Y, Y), chainx)
    Y > 1 && (slicexy .+= kron(matdiag(true, Y, Y, sc=2), chainhopy))

    slicehopxyz = kron(matdiag(true, X*Y, X*Y), hopz)

    grain = kron(matdiag(true, Z, Z), slicexy)
    Z > 1 && (grain .+= kron(matdiag(true, Z, Z, sc=2), slicehopxyz))
    return grain
end

function
main()
    H = weyl(2, 2, 50)
    urange = range(1.00, length=10, stop=30.0)
    hfenergies, hfinter, hfrho = interaction(H, hfcycle; ne=50, urange=urange, T=0.2)
    plot(urange, hfinter)

end

main()
end