module FESetCorotBeamModule

using FinEtools
using ..CrossSectionModule: AbstractCrossSectionType

mutable struct FESetL2CorotBeam{CT} <: AbstractFESet1Manifold{2}
    conn::Array{NTuple{2, FInt}, 1};
    label::FIntVec; 
    crosssection::CT
    A::FFltVec
    I1::FFltVec
    I2::FFltVec
    I3::FFltVec
    J::FFltVec
    x1x2_vector::Vector{FFltVec}
end

function FESetL2CorotBeam(conn::FIntMat, crosssection::CT) where {CT}
    A, J, I1, I2, I3, x1x2_vector = crosssection.parameters(0.0)
    N = size(conn, 1)
    _A = fill(A, N)
    _I1 = fill(I1, N)
    _I2 = fill(I2, N)
    _I3 = fill(I3, N)
    _J = fill(J, N)
    _x1x2_vector = [x1x2_vector for i in 1:N]
    self = FESetL2CorotBeam(NTuple{2, FInt}[], FInt[], crosssection, _A, _I1, _I2, _I3, _J, _x1x2_vector)
    self = fromarray!(self, conn)
    setlabel!(self, 0)
    return self
end

function FESetL2CorotBeam(conn::FIntMat, crosssection::CT, _A, _I1, _I2, _I3, _J, _x1x2_vector) where {CT}
    dummy = FESetL2(conn)
    setlabel!(dummy, 0)
    return FESetL2CorotBeam(dummy.conn, dummy.label, crosssection, _A, _I1, _I2, _I3, _J, _x1x2_vector)
end

function cat(self::T,  other::T) where {T<:FESetL2CorotBeam}
    @assert self.crosssection === other.crosssection "Cannot concatenate sets with distinct cross-sections"
    result = deepcopy(self)
    result.conn = vcat(self.conn, other.conn);
    setlabel!(result, vcat(self.label, other.label))
    result.A = vcat(self.A, other.A);
    result.I1 = vcat(self.I1, other.I1);
    result.I2 = vcat(self.I2, other.I2);
    result.I3 = vcat(self.I3, other.I3);
    result.J = vcat(self.J, other.J);
    result.x1x2_vector = vcat(self.x1x2_vector, other.x1x2_vector);
    return result
end

end # module
