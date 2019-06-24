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

    rho = matdiag(0.20, 2*X, 2*X, sr=1, sc=1, step=2) +
          matdiag(0.80, 2*X, 2*X, sr=2, sc=2, step=2)

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
    X = 41
    H = hubbard(X)

    Trange = 0.00:0.50:1.00
    urange = 0.20:1.20:10.00
#    Trange = 0.1:0.5:4.0
#    urange = 0.25:5.0:12

    rho = matdiag(0.80, 2*X, 2*X, sr=1, sc=1, step=2) +
          matdiag(0.20, 2*X, 2*X, sr=2, sc=2, step=2)

    etotal = zeros(Float64, size(Trange, 1), size(urange, 1))
    neelmag = zeros(Float64, size(Trange, 1), size(urange, 1))
    magnetisation = zeros(Float64, size(Trange, 1), size(urange, 1))

    for (i, T) in enumerate(Trange), (j, u) in enumerate(urange)
        println("Doing T: ", T, " u: ", u)
        etot, dummy, hfrho = hfcycle_new(X, u, T, H, rhostart=rho)

        etotal[i, j] = etot
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
main()
    H = hubbard(50)
    # temperature(H, hcycle; ne=50, u=1.0, Tmin=0.0, Tmax=1.0, runs=100)
    urange = range(0.20, step=0.20, stop=4.0)
    # henergies, hinter, hrho = interaction(H, hcycle; ne=50, urange=urange, T=0.2)
    rho = matdiag(0.49, 100, 100, sc=1, step=2) + matdiag(0.51, 100, 100, sr=2, sc=2, step=2)
   # hfenergies, hfinter, hfrho = interaction(H, hfcycle; rhostart=rho,
   #                                          ne=50, urange=urange, T=0.2)
   #Trange = urange
   hfenergies, hfinter, hfrho = temperature(H, hfcycle; rhostart=rho,
                                            ne=50, u=7.0, Trange=Trange)

    # Make the spin ↓ contribution negative
    hfrho[:, 2:2:end] .*= -1

    display(hfrho[17, :])
    display(sum(hfrho[17, :]))
    plot(urange, sum(hfrho, dims=2))
    xlabel(L"U [t]")
    ylabel(L"\rho_\uparrow - \rho_\downarrow")
#    display(hfrho[2, :])
#    subplot(211)
#    display(sum(hfrho[2, :]))
#    plot(1:50, transpose(reshape(hfrho[2, :], 2, 50)))
#    legend([L"\uparrow", L"\downarrow"])
#
#    subplot(212)
#    display(sum(hfrho[14, :]))
#    plot(1:50, transpose(reshape(hfrho[10, :], 2, 50)))
#    legend([L"\uparrow", L"\downarrow"])
end

end #module