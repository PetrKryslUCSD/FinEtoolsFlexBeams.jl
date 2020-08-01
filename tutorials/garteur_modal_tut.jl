# # Modal analysis of a toy airplane

# ## Description

# This virtual test application is based on the test article  
#     used by the GARTEUR Structures & Materials Action Group 19  
#     which organized a Round Robin exercise where 12 European laboratories  
#     tested a single structure between 1995 and 1997. The benchmark structure   
#     was a laboratory structure built to simulate the dynamic behaviour  
#     of an aeroplane. The structure was initially built for a benchmark  
#     study on experimental modal analysis conducted by the  
#     Structures and Materials Action Group (SM-AG19) of the Group  
#     for Aeronautical Research and Technology in EURope (GARTEUR).  
#         The test-bed was designed and manufactured by ONERA, France.

# ## Goals

# - Show how to construct model from multiple connected beams.

# 

##
# ## Definition of the basic inputs

# The finite element code realize on the basic functionality implemented in this
# package.
using FinEtools

# The material parameters may be defined with the specification of the units.
# The elastic properties are:
E = 70000.0 * phun("MPa") 
nu = 0.31;

# The mass density is
rho = 2700 * phun("kg/m^3")
# This is the characteristic length. The dimensions of the aircraft frame are
# expressed in terms of multiples of this characteristic unit.
L = 0.1*phun("m");


##
# ## Cross-section

# Cross-sectional properties are incorporated in the cross-section property. 
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionRectangle
# Body of the frame.
cs_body = CrossSectionRectangle(s -> 1.5*L, s -> L/2, s -> [1.0, 0.0, 1.0])
# Wing beam.
cs_wing = CrossSectionRectangle(s -> L/10, s -> L, s -> [0.0, 0.0, 1.0])
# Wing drums.
cs_drum = CrossSectionRectangle(s -> L/10, s -> L, s -> [0.0, 0.0, 1.0])
# Vertical part of the tail.
cs_tailv = CrossSectionRectangle(s -> L, s -> L/10, s -> [1.0, 0.0, 1.0])
# Horizontal part of the tail.
cs_tailh = CrossSectionRectangle(s -> L/10, s -> L, s -> [0.0, 0.0, 1.0])
# Massless connectors of the structural parts.
cs_conn = CrossSectionRectangle(s -> L/10, s -> L/10, s -> [1.0, 0.0, 1.0])
# Massless connector between the body and the wings.
cs_connw = CrossSectionRectangle(s -> L/2, s -> L/2, s -> [1.0, 0.0, 1.0])
# Viscoelastic connecting layer.
cs_vconstr = CrossSectionRectangle(s -> L*(1.1/100), s -> L*(76.2/100), s -> [0.0, 0.0, 1.0])

# 
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member
tolerance = L/10000;

# The parts of the mesh can be distinguished based on the label.
meshes = Tuple{FENodeSet, AbstractFESet}[]

# Define the constituent parts of the body of the aircraft.
push!(meshes, frame_member([-9*L 0 0; -8.5*L 0 0], 1, cs_body; label = 1))
push!(meshes, frame_member([-8.5*L 0 0; -8.0*L 0 0], 1, cs_body; label = 1))
push!(meshes, frame_member([-8.0*L 0 0; 0 0 0], 2, cs_body; label = 1))
push!(meshes, frame_member([0 0 0; 6*L 0 0], 2, cs_body; label = 1))

# Define the aluminum parts of the wings.
for i in 1:16
    push!(meshes, frame_member([0 (i-1)*.375*L .805*L;  0 0.375*i*L .805*L], 1, cs_wing; label = 2))
    push!(meshes, frame_member([0 -(i-1)*.375*L .805*L;  0 -0.375*i*L .805*L], 1, cs_wing; label = 2))
end
for i in 1:8
    push!(meshes, frame_member([0 (6*L+(i-1)*.5*L) .805*L; 0 (6*L + 0.5*i*L) .805*L], 1, cs_wing; label = 2))
    push!(meshes, frame_member([0 (-6*L-(i-1)*.5*L) .805*L; 0 (-6*L - 0.5*i*L) .805*L], 1, cs_wing; label = 2))
end

# Define the drums at the ends of the wings.
push!(meshes, frame_member([0 +9.5*L +0.91*L; +2*L +9.5*L +0.91*L], 1, cs_drum; label = 3))
push!(meshes, frame_member([0 +9.5*L +0.91*L; -2*L +9.5*L +0.91*L], 1, cs_drum; label = 3))
push!(meshes, frame_member([0 -9.5*L +0.91*L; +2*L -9.5*L +0.91*L], 1, cs_drum; label = 3))
push!(meshes, frame_member([0 -9.5*L +0.91*L; -2*L -9.5*L +0.91*L], 1, cs_drum; label = 3))

# Define the horizontal and vertical parts of the tail.
push!(meshes, frame_member([-8*L 0 .75*L; -8*L 0 3.35*L], 1, cs_tailv; label = 4))
push!(meshes, frame_member([-8*L 0 3.35*L; -8*L 0 3.75*L], 1, cs_tailv; label = 4))
push!(meshes, frame_member([-8*L 0 3.8*L; -8*L 2*L 3.8*L], 1, cs_tailh; label = 5))
push!(meshes, frame_member([-8*L 0 3.8*L; -8*L -2*L 3.8*L], 1, cs_tailh; label = 5))

# Define the parts of the viscoelastic constraining layer.
for i in 1:16
    push!(meshes, frame_member([-.119*L (i-1)*.375*L .8606*L;  -.119*L 0.375*i*L .8606*L], 1, cs_vconstr; label = 6))
    push!(meshes, frame_member([-.119*L -(i-1)*.375*L .8606*L;  -.119*L -0.375*i*L .8606*L], 1, cs_vconstr; label = 6))
end
for i in 1:5
    push!(meshes, frame_member([-.119*L (6*L+(i-1)*.5*L) .8606*L; -.119*L (6*L+0.5*i*L) .8606*L], 1, cs_vconstr; label = 6))
    push!(meshes, frame_member([-.119*L (-6*L-(i-1)*.5*L) .8606*L; -.119*L (-6*L-0.5*i*L) .8606*L], 1, cs_vconstr; label = 6))
end

# Define the massless connectors between:
# Wing - Wingdrum
push!(meshes, frame_member([0 +9.5*L +0.805*L;0 +9.5*L +0.91*L], 1, cs_conn; label = 7))
push!(meshes, frame_member([0 -9.5*L +0.805*L;0 -9.5*L +0.91*L], 1, cs_conn; label = 7))
# Body-Wing
push!(meshes, frame_member([0 0 0; 0 0 .805*L], 1, cs_connw; label = 7))
# Body-Tail
push!(meshes, frame_member([-8*L 0 0; -8*L 0 .75*L], 1, cs_conn; label = 7))
# Tail-Taildrum
push!(meshes, frame_member([-8*L 0 3.75*L; -8*L 0 3.8*L], 1, cs_conn; label = 7))

# Wing-Viscoelastic Constraining layer
# Middle connector created individually
push!(meshes, frame_member([0 0 .805*L; -.119*L 0 .8606*L], 1, cs_conn; label = 7))

# Connectors created for both wings
for i in 1:16
    push!(meshes, frame_member([0 i*.375*L .805*L; -.119*L 0.375*i*L 0.8606*L], 1, cs_conn; label = 7))
    push!(meshes, frame_member([0 -i*.375*L .805*L; -.119*L -0.375*i*L 0.8606*L], 1, cs_conn; label = 7))
end
for i in 1:6
    push!(meshes, frame_member([0 (6*L+(i-1)*.5*L) .805*L; -.119*L (6*L + 0.5*(i-1)*L) 0.8606*L], 1, cs_conn; label = 7))
    push!(meshes, frame_member([0 (-6*L-(i-1)*.5*L) .805*L; -.119*L (-6*L-0.5*(i-1)*L) 0.8606*L], 1, cs_conn; label = 7))
end

# Massless Sensor Connectors
# Tail Sensors
push!(meshes, frame_member([-8*L 2*L 3.8*L; -(153/20)*L (37/20)*L 3.85*L], 1, cs_conn; label = 8))# 303
push!(meshes, frame_member([-8*L -2*L 3.8*L; -(153/20)*L -(37/20)*L 3.85*L], 1, cs_conn; label = 8))# 301

# Wingdrum Sensors
push!(meshes, frame_member([0 9.5*L .91*L ; 0 9.8*L .96*L], 1, cs_conn; label = 8))# 101
push!(meshes, frame_member([-2*L 9.5*L .91*L ; -1.8*L 9.8*L .96*L], 1, cs_conn; label = 8))# 112
push!(meshes, frame_member([2*L 9.5*L .91*L ; 1.8*L 9.8*L .96*L], 1, cs_conn; label = 8))# 111

push!(meshes, frame_member([0 -9.5*L .91*L ; 0 -9.8*L .96*L], 1, cs_conn; label = 8))# 1
push!(meshes, frame_member([-2*L -9.5*L .91*L ; -1.8*L -9.8*L .96*L], 1, cs_conn; label = 8))# 12
push!(meshes, frame_member([2*L -9.5*L .91*L ; 1.8*L -9.8*L .96*L], 1, cs_conn; label = 8))# 11

# Merge all the meshes of individual parts. This will glued together notes which are in the "same" location.
fens, fesa = mergenmeshes(meshes, tolerance)

# The geometry is visualized in the tutorial garteur_geometry_tut.

##
# ## Material

# Material properties can be now used to create a material: isotropic elasticity model of the `FinEtoolsDeforLinear` package is instantiated.
using FinEtoolsDeforLinear
# The material of the structure is aluminum.
alu = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)
# Material for the massless connectors has the mass density set to zero; otherwise it has the same properties as the aluminum material  of the structure.
massless = MatDeforElastIso(DeforModelRed3D, 0.0, E, nu, 0.0)

# This simple function returns material based on the label of the beam elements.
material(labl) = begin
    if labl == 7 || labl == 8
        return massless
    end
    return alu
end

bungeecoefficient = 4000*phun("N/m");


##
# ## Fields

# Now we start constructing the discrete finite element model.
# We begin by constructing the requisite fields, geometry and displacement.
# These are the so-called "configuration variables", all initialized to 0.
# This is that geometry field.
geom0 = NodalField(fens.xyz)
# This is the displacement field, three unknown displacements per node.
u0 = NodalField(zeros(size(fens.xyz, 1), 3))
# This is the rotation field, three unknown rotations per node are represented
# with a rotation matrix, in total nine numbers. The utility function
# `initial_Rfield`
using FinEtoolsFlexBeams.RotUtilModule: initial_Rfield
Rfield0 = initial_Rfield(fens)
# Here we verify the number of nodes and the number of degrees of freedom in the
# rotation field per node.
@show nents(Rfield0)
@show ndofs(Rfield0)

# Finally, this is the displacement and rotation field for incremental changes,
# incremental displacements and incremental rotations. In total, 6 unknowns per
# node.
dchi = NodalField(zeros(size(fens.xyz, 1), 6))

# There are no support conditions.
applyebc!(dchi)
# The  the number of free
# (unknown) degrees of freedom is equal to the total number of degrees of freedom in the system.
numberdofs!(dchi);

##
# ## Identify support points and locations of sensors

# Suspension points
suspln = selectnode(fens; box = initbox!(Float64[], vec([0.0*L 0.375*L 0.805*L])), inflate = tolerance)
susprn = selectnode(fens; box = initbox!(Float64[], vec([0.0*L -0.375*L 0.805*L])), inflate = tolerance)
suspbn = selectnode(fens; box = initbox!(Float64[], vec([-2.0*L 0.0*L 0.0*L])), inflate = tolerance)

# The sensors at the tip of the left and right wing drum
sensor112n = selectnode(fens; box = initbox!(Float64[], vec([+1.8*L 9.8*L 0.96*L])), inflate = tolerance)
sensor12n = selectnode(fens; box = initbox!(Float64[], vec([+1.8*L -9.8*L .96*L])), inflate = tolerance)
sensor111n = selectnode(fens; box = initbox!(Float64[], vec([-1.8*L 9.8*L 0.96*L])), inflate = tolerance)
sensor11n = selectnode(fens; box = initbox!(Float64[], vec([-1.8*L -9.8*L .96*L])), inflate = tolerance)

# The Connection between the horizontal and vertical tail parts
sensor202n = selectnode(fens; box = initbox!(Float64[], vec([-8*L 0 3.8*L])), inflate = tolerance)

##
# ## Assemble the global discrete system


# For disambiguation we will refer to the stiffness and mass functions by qualifying them with the corotational-beam module, `FEMMCorotBeamModule`.
using FinEtoolsFlexBeams.FEMMCorotBeamModule
CB = FEMMCorotBeamModule

# Thus we can construct the stiffness and mass matrix as follows:
# Note that the finite element machine is the first argument. This provides
# access to the integration domain. The next argument is the geometry field,
# followed by the displacement, rotations, and incremental
# displacement/rotation fields. 
using  SparseArrays

K, M = let
    K = spzeros(dchi.nfreedofs, dchi.nfreedofs)
    M = spzeros(dchi.nfreedofs, dchi.nfreedofs)
    for fes in fesa
        labl  = fes.label[1]
        femm = CB.FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material(labl));
        K += CB.stiffness(femm, geom0, u0, Rfield0, dchi);
        M += CB.mass(femm, geom0, u0, Rfield0, dchi);
    end
    K, M
end

##
# ## Additional concentrated masses.


using LinearAlgebra

using FinEtoolsFlexBeams.FEMMPointMassModule
PM = FEMMPointMassModule

# There is a sensor on the tail.
femmcm1 =  PM.FEMMPointMass(IntegDomain(FESetP1(reshape([sensor202n;], 1, 1)), PointRule()), FFltMat(2*L*L/5*L/5*2*rho*I(3)));

# These are the sensors on the wing drums.
femmcm2 =  PM.FEMMPointMass(IntegDomain(FESetP1(reshape([sensor112n; sensor12n; sensor111n; sensor11n], 4, 1)), PointRule()), FFltMat(0.1*phun("kg")*I(3)));

Mp = PM.mass(femmcm1, geom0, u0, Rfield0, dchi) + PM.mass(femmcm2, geom0, u0, Rfield0, dchi);

##
# ## Bungee supports 


using LinearAlgebra

using FinEtoolsFlexBeams.FEMMPointGroundedSpring
BS = FEMMPointGroundedSpring

# There are three suspension points at the top of the fuselage. We assume that these bungee supports exert only reaction in the vertical direction.
femmbs =  BS.FEMMPointMass(IntegDomain(FESetP1(reshape([suspln; susprn; suspbn;], 3, 1)), PointRule()), 
    FFltMat([bungeecoefficient*[0;0;1]*[0;0;1]' 0*I(3); 0*I(3) 0*I(3)]));

Kb = BS.stiffness(femmbs, geom0, u0, Rfield0, dchi)



# We can compare the size of the stiffness matrix with the number of degrees of
# freedom that are unknown (20).
@show size(K)

##
# ## Solve the free-vibration problem
neigvs = 20
oshift = (2*pi*0.5)^2;

# The Arnoldi algorithm implemented in the well-known `Arpack` package is used
# to solve the generalized eigenvalue problem with the sparse matrices. As is
# common in structural dynamics, we request the smallest eigenvalues in
# absolute value (`:SM`). 
using Arpack
evals, evecs, nconv = eigs(K + oshift * (M + Mp), (M + Mp); nev=neigvs, which=:SM);
# First  we should check that the requested eigenvalues actually converged:
@show nconv == neigvs

# The eigenvalues (i. e. the squares of the angular frequencies) are returned in
# the vector `evals`. The mode shapes constitute the columns of the matrix `evecs`.
@show size(evecs)
# The natural frequencies are obtained from the squares of the angular
# frequencies. We note the use of `sqrt.` which broadcast the square root over
# the array `evals`.
fs = sqrt.([max(0, e - oshift) for e in evals]) / (2 * pi);

##
# ## Comparison of computed and analytical results

# The approximate and analytical frequencies are now reported.
println("Approximate frequencies: $fs [Hz]")



##
# ## Visualize vibration modes

# The animation will show one of the vibration modes overlaid on the undeformed geometry. The configuration during the animation needs to reflect rotations. The function `update_rotation_field!` will update the rotation field given a vibration mode.
using FinEtoolsFlexBeams.RotUtilModule: update_rotation_field!

# The visualization utilities take advantage of the PlotlyJS library.
using PlotlyJS
using FinEtoolsFlexBeams.VisUtilModule: plot_space_box, plot_solid, render, react!, default_layout_3d, save_to_json

# The magnitude of the vibration modes (displacements  and rotations) will be amplified with this scale factor:
scale = 1.5

# This is the mode that will be animated:
mode = 7

# In order to handle variables inside loops correctly, we create a local scope with the `let end` block.
let
    # The extents of the box will be preserved during animation in order to eliminate changes in the viewing parameters.
    tbox = plot_space_box([[-1.2 * radius -1.2 * radius -1.2 * radius]; [+1.2 * radius +1.2 * radius +1.2 * radius]])
    # This is the geometry of the structure without deformation (undeformed). It is displayed as gray, partially transparent.
    tenv0 = plot_solid(fens, fes; x=geom0.values, u=0.0 .* dchi.values[:, 1:3], R=Rfield0.values, facecolor="rgb(125, 155, 125)", opacity=0.3);
    # Initially the plot consists of the box and the undeformed geometry.
    plots = cat(tbox, tenv0; dims=1)
    # Create the layout for the plot. Set the size of the window.
    layout = default_layout_3d(;width=600, height=600, options = Dict(:responsive=>true))
    # Set the aspect mode to get the correct proportions.
    layout[:scene][:aspectmode] = "data"
    # Render the undeformed structure
    pl = render(plots; layout=layout, title="Mode $(mode)")
    sleep(2.115)
    # This is the animation loop. 
    # 1. Distribute a fraction of the selected eigenvector into the incremental displacement/rotation field.
    # 2. Create the deformed configuration by defining displacement field `u1` and rotation field `Rfield1`.
    # 3. Create the plot for the deformed configuration, and add it to the list of plots.
    # 4. Call the `react!` function to update the display. Sleep for a brief period of time to give the display a chance to become current.
    for xscale in scale .* sin.(collect(0:1:89) .* (2 * pi / 21))
        scattersysvec!(dchi, xscale .* evecs[:, mode])
        u1 = deepcopy(u0)
        u1.values .= dchi.values[:, 1:3]
        Rfield1 = deepcopy(Rfield0)
        update_rotation_field!(Rfield1, dchi)
        tenv1 = plot_solid(fens, fes; x=geom0.values, u=dchi.values[:, 1:3], R=Rfield1.values, facecolor="rgb(50, 55, 125)");
        plots = cat(tbox, tenv0, tenv1; dims=1)
        react!(pl, plots, pl.plot.layout)
        sleep(0.115)
    end
    # Save the plot to a Json file. It can be then re-displayed later.
    save_to_json(pl, "deformed_plot.json")
end

# Load the plot from a file.
using FinEtoolsFlexBeams.VisUtilModule: plot_from_json
plot_from_json("deformed_plot.json")
