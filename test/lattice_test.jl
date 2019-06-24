using HFDynamic
function lattice_test()
    ## Graphene
    #tvecs = [sqrt(3) sqrt(3)/2; 0.0 3/2; 0.0 0.0]

    ## Cubic
    tvecs = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    basis = [0.5 0.0; 0.0 0.5; 0.5 0.5]

    function
    peel(pos)
        xs = pos[1]^2
        ys = pos[2]^2
        zs = pos[3]^2

        xsc = (pos[1] + 2)^2
        ysc = (pos[2] + 2)^2
        return xs + ys + zs ≤ 50 && !(xsc + ysc + zs ≤ 50)
    end

    function
    circlehole(pos)
        x = pos[1]
        y = pos[2]
        z = pos[3]

        xs = x^2
        ys = y^2
        zs = z^2

        return xs + ys + zs ≤ 50 && !(-3 ≤ x ≤ 3 && -3 ≤ y ≤ 3 && -10 ≤ z ≤ 10)
    end

    sitecount, lattice = HFDynamic.lattice(tvecs, circlehole; basis=basis, cubesize=10)
    HFDynamic.visualiselattice(lattice; sitecount=sitecount)
end
