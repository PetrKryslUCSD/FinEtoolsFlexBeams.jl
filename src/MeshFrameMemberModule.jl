module MeshFrameMemberModule

# using SimplePCHIP
using Dierckx
using LinearAlgebra: norm
using FinEtools
using ..FinEtoolsFrames.FESetCorotBeamModule: FESetL2CorotBeam

# Mesh of a generally curved beam member given by the location of the
# vertices of the spline curve.
#
# function [fens,fes] = Beam_member(xyz,nL,constructor,Parameters)
#
# xyz = M rows, one vertex location per row;
# nL = Divided into this many elements
# constructor=constructor of the finite element set, default @fe_set_beam3
# Parameters= see fe_set_beam3
# 
function frame_member(xyz, nL, crosssection)
    npts = size(xyz,1);
    s = fill(0.0, npts);
    for j in 2:npts
        s[j] = s[j-1]+norm(xyz[j,:]-xyz[j-1,:]);
    end 
    k = npts < 4 ? npts-1 : 3
    itpx = Spline1D(s, xyz[:, 1]; k = k)
    itpy = Spline1D(s, xyz[:, 2]; k = k)
    itpz = Spline1D(s, xyz[:, 3]; k = k)

    fens, fes = L2block(s[end], nL)

    nxyz = fill(0.0, count(fens), 3)
    for i in 1:count(fens)
        nxyz[i,1] = itpx(fens.xyz[i]);
        nxyz[i,2] = itpy(fens.xyz[i]);
        nxyz[i,3] = itpz(fens.xyz[i]);
    end

    N = count(fes)
    smid = fill(0.0, N)
    c = fes.conn[1]
    smid[1] = norm(fens.xyz[c[2],:]-fens.xyz[c[1],:])/2
    stot = 0.0
    for i in 2:N
        c = fes.conn[i]
        h =  norm(fens.xyz[c[2],:]-fens.xyz[c[1],:])
        smid[i] = smid[i-1] + h/2
        stot = smid[i] + h/2
    end
    
    _A = fill(0.0, N)
    _I1 = fill(0.0, N)
    _I2 = fill(0.0, N)
    _I3 = fill(0.0, N)
    _J = fill(0.0, N)
    _x1x2_vector = [[0.0, 0.0, 1.0] for i in 1:N]
    @show crosssection
    for i in 1:N
        A, J, I1, I2, I3, x1x2_vector = crosssection.parameters(smid[i] / stot)
        _A[i] = A
        _I1[i] = I1
        _I2[i] = I2
        _I3[i] = I3
        _J[i] = J
        _x1x2_vector[i] = deepcopy(x1x2_vector)
    end
    fens.xyz = deepcopy(nxyz)
    fes = FESetL2CorotBeam(connasarray(fes), crosssection, _A, _I1, _I2, _I3, _J, _x1x2_vector)
    return fens, fes
end

end # module
