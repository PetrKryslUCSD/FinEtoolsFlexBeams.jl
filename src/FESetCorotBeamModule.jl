module FESetCorotBeamModule

using FinEtools
using ..CrossSectionModule: AbstractCrossSectionType

mutable struct FESetL2CorotBeam{CT} <: AbstractFESet1Manifold{2}
    conn::Array{NTuple{2, FInt}, 1};
    label::FIntVec; 
    crosssection::CT
end

function FESetL2CorotBeam(conn::FIntMat, crosssection::CT) where {CT}
    self = FESetL2CorotBeam(NTuple{2, FInt}[], FInt[], crosssection)
    self = fromarray!(self, conn)
    setlabel!(self, 0)
    return self
end

end # module
