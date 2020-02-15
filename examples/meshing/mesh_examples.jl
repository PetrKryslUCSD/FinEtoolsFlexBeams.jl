module mesh_examples

using FinEtools
using FinEtoolsFlexBeams.FESetCorotBeamModule: FESetL2CorotBeam
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionCircle, CrossSectionRectangle
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsBeamsVis: plot_nodes, plot_midline, render, plot_space_box

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
    plots = plot_nodes(fens)
    plots = cat(plots, plot_midline(fens, fes; color = "rgb(155, 155, 255)", lwidth = 4), dims = 1)
    # push!(plots, plot_nodes(fens))
    render(plots; aspectratio = [1.0 1.0 4.0])
    true
end # curve_mesh

function argyr_l_frame()
    # Parameters:
    E=71240.0;#MPa
    nu=0.31;# Poisson ratio
    rho=5e-9;
    b=3.0; h=30.0; L=240.0; # cross-sectional dimensions and length of each leg in millimeters

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 1.0, 0.0])

    ## 
    # Choose the mass formulation:
    mass_type=2;

    ## 
    # Reference frequencies
    reffs = [11.2732, 30.5269]
    neigvs = 2;

    # Select the number of elements per half the length.
    xyz = 
    n=8;
    members = []
    push!(members, frame_member([0 0 L; L 0 L], n, cs))
    push!(members, frame_member([L 0 L; L 0 0], n, cs))
    fens, fes = merge_members(members; tolerance = L / 10000);

    plots = cat(plot_space_box([[0 -L/2 0]; [L L/2 L]]),
        plot_nodes(fens), 
        plot_midline(fens, fes; color = "rgb(155, 155, 255)", lwidth = 4), dims = 1)
    render(plots; aspectratio = [1.0 1.0 1.0])
    true
end # argyr_l_frame

function allrun()
    println("#####################################################")
    println("# curve_mesh ")
    curve_mesh()
    println("#####################################################")
    println("# argyr_l_frame ")
    argyr_l_frame()
    return true
end # function allrun

end # module

