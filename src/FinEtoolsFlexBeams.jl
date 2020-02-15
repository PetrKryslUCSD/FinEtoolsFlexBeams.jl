module FinEtoolsFlexBeams

__precompile__(false)

using FinEtools

include("RotUtilModule.jl")
include("CrossSectionModule.jl")
include("FESetCorotBeamModule.jl")
include("MeshFrameMemberModule.jl")
include("FEMMCorotBeamModule.jl")

end # module
