module mesh_examples

using FinEtools
using FinEtoolsFlexBeams.FESetCorotBeamModule: FESetL2CorotBeam
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionCircle
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member
using FinEtoolsBeamsVis: plot_nodes, render

function curve_mesh()
    L = 42
    xyz = [0 0 0;
    0 L/4 L*1/4;
    L/4 L/4 L*2/4;
    L/4 0 L*3/4;
    0 0 L]
    nL = 20
    
    cs = CrossSectionCircle(s -> 5.9910, s -> [0.0, 0.0, 1.0])
    fens, fes = frame_member(xyz, nL, cs)
    plots = [plot_nodes(fens)]
    # push!(plots, plot_nodes(fens))
    render(plots)
    true
end # curve_mesh

function allrun()
    println("#####################################################")
    println("# curve_mesh ")
    curve_mesh()
    return true
end # function allrun

end # module

