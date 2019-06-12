function
lattice(tvecs::Matrix{T}, boundary; basis::Matrix{T}, cubesize::Integer=30) where T <: AbstractFloat
    boxrange = -cubesize:cubesize
    check = []
    for i in boxrange, j in boxrange, k in boxrange
        append!(check, i*tvecs[:, 1] + j*tvecs[:, 2] + k*tvecs[:, 3])
    end
    display(reshape(check, (3, :)))
end