module Wrapper

using LinearAlgebra
using PyPlot
using HFDynamic


function
hubbard(X::Integer)
    t = -1.0
    su = Float64[t 0; 0 t]
    template = matdiag(true, X, X, sc=2)
    template[1, end] = true
    return kron(template, su)
end

function
phasediagram()
    X = 31
    H = hubbard(X)

    Trange = 0.00:0.25:4.0
    urange = 10.00:0.25:10.00
#    Trange = 0.1:0.5:4.0
#    urange = 0.25:5.0:12

    rho = matdiag(0.50, 2*X, 2*X, sr=1, sc=1, step=2) +
          matdiag(0.50, 2*X, 2*X, sr=2, sc=2, step=2)

    neelmag = zeros(Float64, size(Trange, 1), size(urange, 1))
    magnetisation = zeros(Float64, size(Trange, 1), size(urange, 1))
    for (i, T) in enumerate(Trange)
        hfenergies, hfinter, hfrho = interaction(H, hfcycle; rhostart=rho,
                                                 ne=X, urange=urange, T=T)

        display(hfrho)
        # Get spin by subtracting ↑ from ↓ density
        hfrho = hfrho[:, 1:2:end] - hfrho[:, 2:2:end]
        display(hfrho)
        # println(size(magnetisation), " ", size(sum(hfrho, dims=2)))
        magnetisation[i, :] = sum(hfrho, dims=2)
        display(magnetisation)
        # Staggered magnetism to find antiferromagnetic order
        hfrho[:, 1:2:end] .*= -1
        neelmag[i, :] = sum(hfrho, dims=2)
    end
    subplot(211)
    imshow(magnetisation, interpolation=nothing, extent=[urange[1]; urange[end]; Trange[1]; Trange[end]], cmap="gist_ncar")
    xlabel(L"U [t]")
    ylabel(L"\mathrm{T} [t]")
    cbar = colorbar()
    cbar.ax.set_ylabel(L"m")
#    subplot(212)
#    imshow(neelmag, interpolation=nothing, extent=[urange[1]; urange[end]; Trange[1]; Trange[end]], cmap="gist_ncar")
#    xlabel(L"U [t]")
#    ylabel(L"\mathrm{T} [t]")
#    cbar = colorbar()
#    cbar.ax.set_ylabel(L"\mathrm{staggered} m")
end

function
phasediagram_explicit()
    X = 32
    H = hubbard(X)

    Trange = 0.00:0.50:0.00
    urange = 6.50:0.20:6.50
#    Trange = 0.1:0.5:4.0
#    urange = 0.25:5.0:12

    rho = matdiag(0.70, 2*X, 2*X, sr=1, sc=1, step=2) +
          matdiag(0.30, 2*X, 2*X, sr=2, sc=2, step=2)

    etotal = zeros(Float64, size(Trange, 1), size(urange, 1))
    hartree = zeros(Float64, size(Trange, 1), size(urange, 1))
    neelmag = zeros(Float64, size(Trange, 1), size(urange, 1))
    magnetisation = zeros(Float64, size(Trange, 1), size(urange, 1))

    for (i, T) in enumerate(Trange), (j, u) in enumerate(urange)
        println("Doing T: ", T, " u: ", u)
        etot, hfetot, hfrho = hfcycle(X, u, T, H, rhostart=rho)

        etotal[i, j] = etot
        hartree[i, j] = hfetot
        # Get spin by subtracting ↑ from ↓ density
        hfrho = hfrho[1:2:end] - hfrho[2:2:end]
        # println(size(magnetisation), " ", size(sum(hfrho, dims=2)))
        magnetisation[i, j] = sum(hfrho)
        # Staggered magnetism to find antiferromagnetic order
        hfrho[1:2:end] .*= -1
        neelmag[i, j] = sum(hfrho[1:X ÷ 2])
        println("Done")
    end
    plot(urange, transpose(etotal))
#    plot(Trange, hartree)
#    display(magnetisation)
#    display(neelmag)
#    subplot(211)
#    imshow(magnetisation, interpolation=nothing, extent=[urange[1]; urange[end]; Trange[end]; Trange[1]], cmap="gist_ncar")
#    xlabel(L"U [t]")
#    ylabel(L"\mathrm{T} [t]")
#    cbar = colorbar()
#    cbar.ax.set_ylabel(L"m")
#    subplot(212)
#    imshow(neelmag, interpolation=nothing, extent=[urange[1]; urange[end]; Trange[end]; Trange[1]], cmap="gist_ncar")
#    xlabel(L"U [t]")
#    ylabel(L"\mathrm{T} [t]")
#    cbar = colorbar()
#    cbar.ax.set_ylabel(L"\mathrm{staggered } m")
end

function
teststart()
    X = 32
    H = hubbard(X)

    T = 0.00
    u = 4.50

    startrange = 0:0.01:0.5
    etotal = zeros(Float64, size(startrange, 1))
    hartree = zeros(Float64, size(startrange, 1))
    neelmag = zeros(Float64, size(startrange, 1))
    magnetisation = zeros(Float64, size(startrange, 1))

    for (i, s) in enumerate(startrange)
        rho = matdiag(1-s, 2*X, 2*X, sr=1, sc=1, step=2) +
              matdiag(s, 2*X, 2*X, sr=2, sc=2, step=2)


        etot, hfetot, hfrho = hfcycle(X, u, T, H, rhostart=rho)

        etotal[i] = etot
        hartree[i] = hfetot
        # Get spin by subtracting ↑ from ↓ density
        hfrho = hfrho[1:2:end] - hfrho[2:2:end]
        # println(size(magnetisation), " ", size(sum(hfrho, dims=2)))
        magnetisation[i] = sum(hfrho)
        # Staggered magnetism to find antiferromagnetic order
        hfrho[1:2:end] .*= -1
        neelmag[i] = sum(hfrho[1:X ÷ 2])
        println("Done")
    end
    plot(startrange, etotal)
end

end #module