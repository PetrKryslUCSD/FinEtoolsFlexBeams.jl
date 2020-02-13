module FinEtoolsFrames

using FinEtools

mutable struct FESetL2CorotBeam <: AbstractFESet1Manifold{2}
    conn::Array{NTuple{2, FInt}, 1};
    label::FIntVec; 
    function FESetL2CorotBeam(conn::FIntMat)
        self = new(NTuple{2, FInt}[], FInt[])
        self = fromarray!(self, conn)
        setlabel!(self, 0)
        return self
    end
end

include("MeshFrameMemberModule.jl")

end # module
