module FEMMCorotBeamModule

using LinearAlgebra: norm, Transpose, mul!
using FinEtools
using FinEtools.IntegDomainModule: IntegDomain
import FinEtoolsDeforLinear.MatDeforElastIsoModule: MatDeforElastIso
using ..FESetCorotBeamModule: FESetL2CorotBeam, local_frame_and_def!, local_mass!, local_stiffness!, natural_forces!, local_geometric_stiffness!, MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA

"""
    FEMMDeforLinear{MR<:AbstractDeforModelRed,  S<:AbstractFESet, F<:Function, M<:AbstractMatDeforLinearElastic} <: FEMMCorotBeam

Class for linear deformation finite element modeling machine.
"""
mutable struct FEMMCorotBeam{S<:AbstractFESet, F<:Function} <: AbstractFEMM
    integdomain::IntegDomain{S, F} # integration domain data
    material::MatDeforElastIso # material object
    _ecoords0::FFltMat
    _ecoords1::FFltMat
    _edisp1::FFltMat
    _dofnums::FIntMat
    _F0::FFltMat
    _Ft::FFltMat
    _FtI::FFltMat
    _FtJ::FFltMat
    _Te::FFltMat
    _elmat::FFltMat
    _elmatTe::FFltMat
    _aN::FFltMat
    _dN::FFltVec
    _DN::FFltMat
    _PN::FFltVec
    _RI::FFltMat
    _RJ::FFltMat
end

function _buffers(self)
    return self._ecoords0, self._ecoords1, self._edisp1, self._dofnums, self._F0, self._Ft, self._FtI, self._FtJ, self._Te, self._elmat, self._elmatTe, self._aN, self._dN, self._DN, self._PN, self._RI, self._RJ
end

function _transfmat!(Te, Ft)
    Te[1:3, 1:3] = Te[4:6, 4:6] = Te[7:9, 7:9] = Te[10:12, 10:12] = Ft
    return Te
end

function FEMMCorotBeam(integdomain::IntegDomain{S, F}, material::MatDeforElastIso) where {S<:FESetL2CorotBeam, F<:Function}
    _ecoords0 = fill(0.0, 2, 3)
    _ecoords1 = fill(0.0, 2, 3)
    _edisp1 = fill(0.0, 2, 3)
    _dofnums = zeros(FInt, 1, 12); 
    _F0 = fill(0.0, 3, 3)
    _Ft = fill(0.0, 3, 3)
    _FtI = fill(0.0, 3, 3)
    _FtJ = fill(0.0, 3, 3)
    _Te = fill(0.0, 12, 12)
    _elmat = fill(0.0, 12, 12)
    _elmatTe = fill(0.0, 12, 12)
    _aN = fill(0.0, 6, 12)
    _dN = fill(0.0, 6)
    _DN = fill(0.0, 6, 6)
    _PN = fill(0.0, 6)
    _RI = fill(0.0, 3, 3)
    _RJ = fill(0.0, 3, 3)
    return FEMMCorotBeam(integdomain, material, _ecoords0, _ecoords1, _edisp1, _dofnums, _F0, _Ft, _FtI, _FtJ, _Te, _elmat, _elmatTe, _aN, _dN, _DN, _PN, _RI, _RJ)
end


"""
    mass(self::FEMMCorotBeam,  assembler::A,
      geom::NodalField{FFlt},
      u::NodalField{T}) where {A<:AbstractSysmatAssembler, T<:Number}

Compute the consistent mass matrix

This is a general routine for the abstract linear-deformation  FEMM.
"""
function mass(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}; mass_type=MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA) where {ASS<:AbstractSysmatAssembler, T<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums, F0, Ft, FtI, FtJ, Te, elmat, elmatTe, aN, dN, DN, PN, R1I, R1J = _buffers(self)
    rho = massdensity(self.material)
    A, I1, I2, I3, x1x2_vector = fes.A, fes.I1, fes.I2, fes.I3, fes.x1x2_vector
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), dchi.nfreedofs, dchi.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        R1I[:] .= Rfield1.values[fes.conn[i][1], :];
        R1J[:] .= Rfield1.values[fes.conn[i][2], :];
        fill!(elmat,  0.0); # Initialize element matrix
        L1, Ft, dN = local_frame_and_def!(Ft, dN, F0, FtI, FtJ, ecoords0, x1x2_vector[i], ecoords1, R1I, R1J);
        L0 = norm(ecoords0[2,:]-ecoords0[1,:]); 
        local_mass!(elmat, A[i], I1[i], I2[i], I3[i], rho, L0, mass_type);
        _transfmat!(Te, Ft)
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function mass(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}; mass_type=MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return mass(self, assembler, geom0, u1, Rfield1, dchi; mass_type = mass_type);
end

"""
    stiffness(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the material stiffness matrix.
"""
function stiffness(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums, F0, Ft, FtI, FtJ, Te, elmat, elmatTe, aN, dN, DN, PN, R1I, R1J = _buffers(self)
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)
    A, I2, I3, J, x1x2_vector = fes.A, fes.I2, fes.I3, fes.J, fes.x1x2_vector
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), dchi.nfreedofs, dchi.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        R1I[:] .= Rfield1.values[fes.conn[i][1], :];
        R1J[:] .= Rfield1.values[fes.conn[i][2], :];
        fill!(elmat,  0.0); # Initialize element matrix
        L1, Ft, dN = local_frame_and_def!(Ft, dN, F0, FtI, FtJ, ecoords0, x1x2_vector[i], ecoords1, R1I, R1J);
        local_stiffness!(elmat, E, G, A[i], I2[i], I3[i], J[i], L1, aN, DN);
        _transfmat!(Te, Ft)
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function stiffness(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return stiffness(self, assembler, geom0, u1, Rfield1, dchi);
end


"""
    geostiffness(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}

Compute the geometric stiffness matrix.
"""
function geostiffness(self::FEMMCorotBeam, assembler::ASS, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {ASS<:AbstractSysmatAssembler, T<:Number}
    fes = self.integdomain.fes
    ecoords0, ecoords1, edisp1, dofnums, F0, Ft, FtI, FtJ, Te, elmat, elmatTe, aN, dN, DN, PN, R1I, R1J = _buffers(self)
    E = self.material.E
    G = E / 2 / (1 + self.material.nu)
    A, I2, I3, J, x1x2_vector = fes.A, fes.I2, fes.I3, fes.J, fes.x1x2_vector
    startassembly!(assembler, size(elmat, 1), size(elmat, 2), count(fes), dchi.nfreedofs, dchi.nfreedofs);
    for i = 1:count(fes) # Loop over elements
        gathervalues_asmat!(geom0, ecoords0, fes.conn[i]);
        gathervalues_asmat!(u1, edisp1, fes.conn[i]);
        ecoords1 .= ecoords0 .+ edisp1
        R1I[:] .= Rfield1.values[fes.conn[i][1], :];
        R1J[:] .= Rfield1.values[fes.conn[i][2], :];
        fill!(elmat,  0.0); # Initialize element matrix
        L1, Ft, dN = local_frame_and_def!(Ft, dN, F0, FtI, FtJ, ecoords0, x1x2_vector[i], ecoords1, R1I, R1J);
        natural_forces!(PN, E, G, A[i], I2[i], I3[i], J[i], L1, dN, DN)
        local_geometric_stiffness!(elmat, A[i], I2[i], I3[i], PN, L1)
        _transfmat!(Te, Ft)
        mul!(elmatTe, elmat, Transpose(Te))
        mul!(elmat, Te, elmatTe)
        gatherdofnums!(dchi, dofnums, fes.conn[i]); # degrees of freedom
        assemble!(assembler, elmat, dofnums, dofnums); 
    end # Loop over elements
    return makematrix!(assembler);
end

function geostiffness(self::FEMMCorotBeam, geom0::NodalField{FFlt}, u1::NodalField{T}, Rfield1::NodalField{T}, dchi::NodalField{T}) where {T<:Number}
    assembler = SysmatAssemblerSparseSymm();
    return geostiffness(self, assembler, geom0, u1, Rfield1, dchi);
end

end # module
