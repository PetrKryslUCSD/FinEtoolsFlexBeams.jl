# # Modal analysis of Argyris frame: effect of prestress

# ## Description

# Vibration analysis of a L-shaped frame under a loading. 
# The fundamental vibration frequency depends on the prestress force.

# ## Goals

# - 

# The finite element code relies on the basic functionality implemented in this
# package.
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexBeams.RotUtilModule:  update_rotation_field!
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexBeams.RotUtilModule: initial_Rfield
using FinEtoolsFlexBeams.FEMMCorotBeamModule: FEMMCorotBeam
using FinEtoolsFlexBeams.FEMMCorotBeamModule
CB = FEMMCorotBeamModule

# Parameters:
E = 71240.0 * phun("MPa")
nu = 0.31; # Poisson ratio
rho = 5000 * phun("kg/m^3");
# cross-sectional dimensions and length of each leg in millimeters
b = 0.6 * phun("mm"); h = 30.0 * phun("mm"); L = 240.0 * phun("mm"); 
# Magnitude of the total applied force, Newton
magn = 1e-5 * phun("N");

# Cross-sectional properties
cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 1.0, 0.0])

# Select the number of elements per leg.
n = 8;
members = Tuple{FENodeSet, AbstractFESet}[]
push!(members, frame_member([0 0 L; L 0 L], n, cs))
push!(members, frame_member([L 0 L; L 0 0], n, cs))
fens, fes = merge_members(members; tolerance = L / 10000)

# Material properties
material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

# Construct the requisite fields, geometry and displacement
# Initialize configuration variables
geom0 = NodalField(fens.xyz)
u0 = NodalField(zeros(size(fens.xyz,1), 3))
Rfield0 = initial_Rfield(fens)
dchi = NodalField(zeros(size(fens.xyz,1), 6))

# Apply EBC's
l1 = selectnode(fens; box = [0 0 0 0 L L], tolerance = L / 10000)
for i in [1, 2, 3, 4, 5, 6]
    setebc!(dchi, l1, true, i)
end
applyebc!(dchi)
numberdofs!(dchi);

# Assemble the global discrete system
femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
K = CB.stiffness(femm, geom0, u0, Rfield0, dchi);
M = CB.mass(femm, geom0, u0, Rfield0, dchi);

tipn = selectnode(fens; box=[L L 0 0  0 0], tolerance=L/n/1000)[1]
loadbdry = FESetP1(reshape([tipn], 1, 1))
lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
fi = ForceIntensity(FFlt[-magn, 0, 0, 0, 0, 0]);

F = CB.distribloads(lfemm, geom0, dchi, fi, 3);

# Solve for the displacement under @show the static load
scattersysvec!(dchi, K\F);

# Update deflections so that the initial stress can be computed. First the displacements:
u1 = deepcopy(u0)
u1.values .= dchi.values[:, 1:3]
# Then the rotations:
Rfield1 = deepcopy(Rfield0)
update_rotation_field!(Rfield1, dchi)

# The static deflection is now used to compute the internal forces
# which in turn lead to the geometric stiffness.
Kg = CB.geostiffness(femm, geom0, u1, Rfield1, dchi);

using Arpack
neigvs = 4

lfp = linearspace(0.0, 68000.0, 400)
fsp = let
    fsp = []
    for load_factor in lfp
        evals, evecs, nconv = eigs(K + load_factor .* Kg, M; nev=neigvs, which=:SM);
        f = evals[1] > 0 ? sqrt(evals[1]) / (2 * pi) : 0;
        push!(fsp, f)
    end
    fsp
end

lfm = linearspace(-109000.0, 0.0, 400)
fsm = let
    fsm = []
    for load_factor in lfm
        evals, evecs, nconv = eigs(K + load_factor .* Kg, M; nev=neigvs, which=:SM);
        f = evals[1] > 0 ? sqrt(evals[1]) / (2 * pi) : 0;
        push!(fsm, f)
    end
    fsm
end

using PlotlyJS

tcp = scatter(;x=cat(lfp, lfm; dims=1), y=cat(fsp, fsm; dims=1), mode="markers", name = "Fundamental frequency", line_color = "rgb(15, 15, 15)")
plots = cat(tcp; dims = 1)
layout = Layout(;width=500, height=500, xaxis=attr(title="Loading factor P", zeroline=true), yaxis=attr(title="Frequency(P) [Hz]", zeroline=true))
pl = plot(plots, layout)
display(pl)

##
# ## Visualize some fundamental mode shapes

using FinEtoolsFlexBeams.VisUtilModule: plot_space_box, plot_solid, render, react!, default_layout_3d, save_to_json
scale = 0.005

vis(loading_factor, evec) = let
    tbox = plot_space_box(reshape(inflatebox!(boundingbox(fens.xyz), 0.3 * L), 2, 3))
    tenv0 = plot_solid(fens, fes; x=geom0.values, u=0.0 .* dchi.values[:, 1:3], R=Rfield0.values, facecolor="rgb(125, 155, 125)", opacity=0.3);
    plots = cat(tbox, tenv0; dims=1)
    layout = default_layout_3d(;width=600, height=600)
    layout[:scene][:aspectmode] = "data"
    pl = render(plots; layout=layout, title = "Loading factor $(loading_factor)")
    for xscale in scale .* sin.(collect(0:1:89) .* (2 * pi / 21))
        scattersysvec!(dchi, xscale .* evec)
        u1 = deepcopy(u0)
        u1.values .= dchi.values[:, 1:3]
        Rfield1 = deepcopy(Rfield0)
        update_rotation_field!(Rfield1, dchi)
        tenv1 = plot_solid(fens, fes; x=geom0.values, u=dchi.values[:, 1:3], R=Rfield1.values, facecolor="rgb(50, 55, 125)");
        plots = cat(tbox, tenv0, tenv1; dims=1)
        react!(pl, plots, pl.plot.layout)
        sleep(0.115)
    end
end

loading_factor = 60000
evals, evecs, nconv = eigs(K + loading_factor .* Kg, M; nev=neigvs, which=:SM);
vis(loading_factor, evecs[:, 1])

loading_factor = -50000
evals, evecs, nconv = eigs(K + loading_factor .* Kg, M; nev=neigvs, which=:SM);
vis(loading_factor, evecs[:, 1])

loading_factor = -100000
evals, evecs, nconv = eigs(K + loading_factor .* Kg, M; nev=neigvs, which=:SM);
vis(loading_factor, evecs[:, 1])

true
