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
using PlotlyJS
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
    nodes = scatter3d(;x=fens.xyz[:, 1],y=fens.xyz[:, 2], z=fens.xyz[:, 3], mode="lines",
        marker=attr(color="#1f77b4", size=12, symbol="circle",
            line=attr(color="rgb(0,0,0)", width=3)),
        line=attr(color="#1f77b4", width=4))
    layout = Layout(autosize=false, width=500, height=500,
        margin=attr(l=0, r=0, b=0, t=65))
    display(plot([nodes], layout))
end
end
using .m2
m2.test()
