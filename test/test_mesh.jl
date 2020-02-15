module m1
using FinEtools
using FinEtoolsFlexBeams.FESetCorotBeamModule: FESetL2CorotBeam
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionCircle
using Test
function test()
    cs = CrossSectionCircle(s -> 5.9910, s -> [0.0, 0.0, 1.0])
    fes = FESetL2CorotBeam([1 2; 2 3; 3 4], cs)
end
end
using .m1
m1.test()

module m2
using FinEtools
using FinEtoolsFlexBeams.FESetCorotBeamModule: FESetL2CorotBeam
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionCircle
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member
using Test
function test()
    L = 42
    xyz = [0 0 0;
    0 L/4 L*1/4;
    L/4 L/4 L*2/4;
    L/4 0 L*3/4;
    0 0 L]
    nL = 20
    
    cs = CrossSectionCircle(s -> 5.9910, s -> [0.0, 0.0, 1.0])
    fens, fes = frame_member(xyz, nL, cs)
end
end
using .m2
m2.test()
