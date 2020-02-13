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
    # itpx = interpolate(s, xyz[:, 1])
    # itpy = interpolate(s, xyz[:, 2])
    # itpz = interpolate(s, xyz[:, 3])
    itpx = Spline1D(s, xyz[:, 1])
    itpy = Spline1D(s, xyz[:, 2])
    itpz = Spline1D(s, xyz[:, 3])

    fens, fes = L2block(s[end], nL)

    nxyz = fill(0.0, count(fens), 3)
    for i in 1:count(fens)
        nxyz[i,1] = itpx(fens.xyz[i]);
        nxyz[i,2] = itpy(fens.xyz[i]);
        nxyz[i,3] = itpz(fens.xyz[i]);
    end
  
    fens.xyz = deepcopy(nxyz)
    fes = FESetL2CorotBeam(connasarray(fes), crosssection)
    return fens, fes
end

end # module
