"""
    matdiag(diag, nr, nc; <keyword arguments>)

Create Matrix with number `vdiag` on the super- or subdiagonals and `vndiag` 
in the rest.

# Arguments
- `diag::Number`: `Number` to write into created super- or subdiagonal
- `nr::Integer`: Number of rows
- `nc::Integer`: Number of columns
- `sr::Integer=1`: Starting row of blocks
- `sc::Integer=1`: Starting column of blocks
- `er::Integer=nblocks(M, 1)`: Ending row of blocks
- `ec::Integer=nblocks(M, 2)`: Ending column of blocks
- `step::Integer=1`: Stepwidth on the diagonal of blocks
- `ndiag::Number`: `Number` to write into everything but the specified diagonal

# Examples
```jldoctest
julia> matdiag(true, 5, 5, sr=2, ec=3)
5×5 Array{Bool,2}:
 false  false  false  false  false
  true  false  false  false  false
 false   true  false  false  false
 false  false   true  false  false
 false  false  false  false  false
```
"""
function
matdiag(diag::Number, nr::Integer, nc::Integer; 
        sr::Integer=1, sc::Integer=1,
        er::Integer=nr, ec::Integer=nc,
        step::Integer=1, ndiag::Number=zero(diag))

    M = fill(ndiag, nr, nc)
    matdiag!(diag, M, sr=sr, sc=sc, er=er, ec=ec, step=step)
    return M
end


"""
    matdiag!(diag, M; <keyword arguments>)

Put Number `v` on the super- or subdiagonals of Matrix `M`.

# Arguments
- `diag`: `Number` to write into `M`s diagonal
- `M::Matrix`: Matrix to be written to
- `sr::Integer=1`: Starting row of blocks
- `sc::Integer=1`: Starting column of blocks
- `er::Integer=nblocks(M, 1)`: Ending row of blocks
- `ec::Integer=nblocks(M, 2)`: Ending column of blocks
- `step::Integer=1`: Stepwidth on the diagonal of blocks

# Examples
```jldoctest
julia> A = Matrix(fill(false, 6, 6))
6×6 Array{Bool,2}:
 false  false  false  false  false  false
 false  false  false  false  false  false
 false  false  false  false  false  false
 false  false  false  false  false  false
 false  false  false  false  false  false
 false  false  false  false  false  false

julia> matdiag!(A, s, sc=2, step=2)
julia> A
6x6 Array{Bool,2}:
 false   true  false  false  false  false
 false  false  false  false  false  false
 false  false  false   true  false  false
 false  false  false  false  false  false
 false  false  false  false  false   true
 false  false  false  false  false  false
```
"""
function
matdiag!(diag::Number, M::Matrix;
         sr::Integer=1, sc::Integer=1,
         er::Integer=size(M, 1), ec::Integer=size(M, 2),
         step::Integer=1)

    (sr > er || sc > ec) &&
        error("Starting row or starting column too large. ",
              "sr=", sr, "<", er, " sc=", sc, "<", ec)
    
    while sr ≤ er && sc ≤ ec
        M[sr, sc] = diag
        sr += step
        sc += step
    end
    return nothing
end