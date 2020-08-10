module MeshFrameMemberModule

# using SimplePCHIP
using Dierckx
using LinearAlgebra: norm
using FinEtools
using ..FinEtoolsFlexBeams.FESetCorotBeamModule: FESetL2CorotBeam
using ..FinEtoolsFlexBeams.FESetCorotBeamModule: cat

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
function frame_member(xyz, nL, crosssection; label = 0)
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
    fens.xyz = deepcopy(nxyz)
    
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
    _A2s = fill(0.0, N)
    _A3s = fill(0.0, N)
    _x1x2_vector = [[0.0, 0.0, 1.0] for i in 1:N]
    _dimensions = [[0.0] for i in 1:N]
    for i in 1:N
        par = crosssection.parameters(smid[i] / stot)
        _A[i] = par.A
        _I1[i] = par.I1
        _I2[i] = par.I2
        _I3[i] = par.I3
        _J[i] = par.J
        _A2s[i] = par.A2s
        _A3s[i] = par.A3s
        _x1x2_vector[i] = deepcopy(par.x1x2_vector)
        _dimensions[i] = deepcopy(par.dimensions)
    end
    fes = FESetL2CorotBeam(connasarray(fes), crosssection, _A, _I1, _I2, _I3, _J, _A2s, _A3s, _x1x2_vector, _dimensions)
    setlabel!(fes, label)
    return fens, fes
end

function fuse_members(members; tolerance = 0.001)
    Meshes = Array{Tuple{FENodeSet, AbstractFESet},1}()
    for m in members
        push!(Meshes, m)
    end
    return mergenmeshes(Meshes, tolerance);
end

function merge_members(members; tolerance = 0.001)
    fens, allfes = fuse_members(members; tolerance = tolerance)
    fes = allfes[1]
    for ofes in allfes[2:end]
        fes = cat(fes, ofes)
    end
    return fens, fes
end

end # module
