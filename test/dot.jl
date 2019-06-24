module Wrapper
using BenchmarkTools
using LinearAlgebra

include("tools.jl")
include("hartree.jl")
include("eval.jl")

using .Tools
using .Hartree
using .Eval

function
quantumdot(X::Integer, dotpos::Integer)
    t = 1.0
    su = Float64[t 0; 0 t]
    chainx = kron(matdiag(true, X, X, sc=2), su)
    up = 2*dotpos-1
    down = 2*dotpos

    on = 0.1t    
    v = 0.2t
    chainx[up, up] = on
    chainx[down, down] = on
    chainx[up - 2, up] = v
    chainx[down - 2, down] = v
    return chainx
end

function
main()
    X = 50
    dotpos = X
    H = quantumdot(X, dotpos)
    #temperature(H, dotcycle, dotpos; ne=50, u=1.0, Tmin=0.0, Tmax=1.0, runs=100)
    interaction(H, dotcycle, dotpos; ne=50, umin=0.01, umax=2.0, T=0.2, runs=100)
    #electrons(H, dotcycle, dotpos; nemin=10, nemax=80, u=1.0, T=0.2, nestep=1)
end
end
