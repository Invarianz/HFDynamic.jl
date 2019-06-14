"""
    lattice(tvecs, boundary; basis, cubesize)

Create 1d, 2d or 3d lattice from translational vectors and an additional set of
optional basis vectors. The returned lattice is ordered by basis sites. E.g. if
the lattice has 10 sites 5 belonging to basis type A and 5 to basis type B the
`lattice[:, 1:5]` are type A sites and `lattice[:, 6:10]` type B sites. Correspondingly
`sitecount` reads `[5 5]`.

# Arguments
- `tvecs::Matrix{<:AbstractFloat}`: translational vectors in column major order
- `boundary`: function determining the lattice boundary, needs to take a position vector with 3 components and return true or false
- `basis::Matrix{<:AbstractFloat}`: *optional* basis vectors in column major order
- `cubesize::Int64`: *optional* Size of the chain (1d), rectangle (2d) or cube (3d) cut with the boundary function
"""
function
lattice(tvecs::Matrix{T}, boundary; basis::Matrix{T}=zeros(T, 3, 1), 
        cubesize::Int64=30) where T <: AbstractFloat

    # Add the zero basis vector to the problem if a basis was defined
    if basis != zeros(T, 3, 1)
        basis = hcat(zeros(T, 3, 1), basis)
    end
    # The complete lattice with basis atoms
    carvedlattice = zeros(T, 0)

    # How many sites belong to which basis vector, holds only one number
    # if no basis is specified (number of sites)
    sitecount = zeros(Int64, size(basis, 2))

    if size(tvecs, 2) == 3
        # 3d system (3 translational vectors)
        for (i, bpos) in enumerate(eachcol(basis))
            sitecount[i] = _cube!(carvedlattice, bpos, tvecs, boundary, cubesize)
        end
    elseif size(tvecs, 2) == 2
        # 2d system (2 translational vectors)
        for (i, bpos) in enumerate(eachcol(basis))
            sitecount[i] = _film!(carvedlattice, bpos, tvecs, boundary, cubesize)
        end
    elseif size(tvecs, 1) == 1
        # 1d system (1 translational vector)
        for (i, bpos) in enumerate(eachcol(basis))
            sitecount[i] = _chain!(carvedlattice, bpos, tvecs, boundary, cubesize)
        end
    end

    return sitecount, reshape(carvedlattice, (3, :))
end

function
_cube!(carvedlattice, bpos, tvecs::Matrix{T}, boundary, 
       cubesize::Integer) where T <: AbstractFloat
    boxrange = -cubesize:cubesize
    sitecount = 0
    for i in boxrange, j in boxrange, k in boxrange
        pos = i*tvecs[:, 1] + j*tvecs[:, 2] + k*tvecs[:, 3] + bpos
        if boundary(pos)
            append!(carvedlattice, pos)
            sitecount += 1
        end
    end
    return sitecount
end

function
_film!(carvedlattice, bpos, tvecs::Matrix{T}, boundary,
       cubesize::Integer) where T <: AbstractFloat
    boxrange = -cubesize:cubesize
    sitecount = 0
    for i in boxrange, j in boxrange, k in boxrange
        pos = i*tvecs[:, 1] + j*tvecs[:, 2] + bpos
        if boundary(pos)
            append!(carvedlattice, pos)
            sitecount += 1
        end
    end
    return sitecount
end

function
_chain!(carvedlattice, bpos, tvecs::Matrix{T}, boundary,
        cubesize::Integer) where T <: AbstractFloat
    boxrange = -cubesize:cubesize
    sitecount = 0
    for i in boxrange, j in boxrange, k in boxrange
        pos = i*tvecs[:, 1] + bpos
        if boundary(pos)
            append!(carvedlattice, pos)
            sitecount += 1
        end
    end
    return sitecount
end

"""
    visualiselattice(lattice; sitecount)

Visualise the lattice in a plot. Lattice needs to be ordered column wise i.e.
(3 x #sites). `sitecount` is an array listing the number of sites belonging to the
corresponding basis. I.e. if lattice has 10 sites whereas 5 belong to basis type A
and 5 to basis type B, sitecount should read `[5 5]`.
"""
function
visualiselattice(lattice::Matrix{<:AbstractFloat};
                 sitecount::Array{Int64, 1}=[size(lattice, 2)])

    plotlattice = Array{Point{3, Float32},1}(undef, size(lattice, 2))
    for (i, pos) in enumerate(eachcol(lattice))
        plotlattice[i] = Point3f0(pos)
    end

    color = Array{Symbol, 1}(undef, size(lattice, 2))
    offset = 1
    colors = [:red, :blue, :green, :yellow, :brown]
    for (i, sites) in enumerate(sitecount)
        color[offset:(offset-1)+sites] .= colors[i]
        offset += sites
    end
    meshscatter(plotlattice, color=color)
end