# # GARTEUR SM-AG19 Testbed: Harmonic Vibration Analysis 

# ## Description

# This virtual test application is based on the test article  
# used by the GARTEUR Structures & Materials Action Group 19  
# which organized a Round Robin exercise where 12 European laboratories  
# tested a single structure between 1995 and 1997. The benchmark structure   
# was a laboratory structure built to simulate the dynamic behaviour  
# of an aeroplane. The structure was initially built for a benchmark  
# study on experimental modal analysis conducted by the  
# Structures and Materials Action Group (SM-AG19) of the Group  
# for Aeronautical Research and Technology in EURope (GARTEUR).  
# The test-bed was designed and manufactured by ONERA, France.

# ![](IMAC97photo.jpg)

# ### References

# [1] Ground Vibration Test Techniques, compiled by A Gravelle, GARTEUR
# Structures & Materials Action Group 19 Technical report TP-115, 1999.
# [2] Etienne Balmes, Jan R. Wright, GARTEUR GROUP ON GROUND VIBRATION TESTING |
# RESULTS FROM THE TEST OF A SINGLE STRUCTURE BY 12 LABORATORIES IN EUROPE,
# Proceedings of DETC'97, 1997 ASME Design Engineering Technical Conferences,
# September 14-17, 1997, Sacramento, California.
# [3] 3M(TM) Viscoelastic Damping Polymer 112 Series,  Technical Data, May 2017.

# ## Goals

# - Show how to construct model from multiple connected beams.
# - Demonstrate the use of massless connectors.
# - Demonstrate the use of point masses.
# - Demonstrate the use of grounded springs.  
# - Illustrate verification of the solution of the free vibration problem. 

##
# ## Geometry of the testbed airplane.

# It was a rather simple structure which was reasonably dynamically
# representative of a simple airplane structure. It was composed of several beams
# simulating a fuselage with wings and a tail. Wing tip drums allowed to adjust
# bending and torsion frequencies similarly to airplane ones, with some very
# close modal frequencies. 

# ![](garteur-geom.png)

# The beam finite element code relies on the basic functionality implemented in this
# package.
using FinEtools

# This is the characteristic length. The dimensions of the aircraft frame are
# expressed in terms of multiples of this characteristic unit.
L = 0.1*phun("m");

##
# ## Cross-section

# Cross-sectional properties are incorporated in the cross-section property.
# There are several rectangular cross-sections in the model: the fuselage, the
# wing, the tail. There are also three massless connectors: connections between
# the fuselage and the wing, between the wing structure and the viscoelastic
# damping layer, and between the fuselage and the tail.
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
cs_conn = CrossSectionRectangle(s -> L/5, s -> L/5, s -> [1.0, 0.0, 1.0])
# Massless connector between the body and the wings.
cs_connw = CrossSectionRectangle(s -> L/2, s -> L, s -> [0.0, 1.0, 0.0])
# Massless connector between the body and the tail.
cs_connt = CrossSectionRectangle(s -> L, s -> L/5, s -> [1.0, 0.0, 1.0])
# Viscoelastic connecting layer.
cs_vconstr = CrossSectionRectangle(s -> L*(1.1/100), s -> L*(76.2/100), s -> [0.0, 0.0, 1.0])

# ## Mesh 

# We shall use this utility function to generate the mesh of the individual
# parts. This will result in a number of separate meshes for the members. These
# separate meshes will then be glued together (merged) based on the tolerance on
# the location of the nodes.
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member
tolerance = L/10000;
# Number of intervals from 0.25*L to 8.5*L (the extent of the constraining plate).
nc = 8
meshes = Tuple{FENodeSet, AbstractFESet}[]

# Define the constituent parts of the body of the aircraft.
push!(meshes, frame_member([-9*L 0 0; -8.5*L 0 0], 1, cs_body; label = 1))
push!(meshes, frame_member([-8.5*L 0 0; -8.0*L 0 0], 1, cs_body; label = 1))
push!(meshes, frame_member([-8.0*L 0 0; -2.0*L 0 0], 2, cs_body; label = 1))
push!(meshes, frame_member([-2.0*L 0 0; 0 0 0], 1, cs_body; label = 1))
push!(meshes, frame_member([0 0 0; 6*L 0 0], 2, cs_body; label = 1))

# Define the aluminum parts of the wings.
push!(meshes, frame_member([0 0 0.805*L;  0 0.25*L 0.805*L], 1, cs_wing; label = 2))
push!(meshes, frame_member([0 0 0.805*L;  0 -0.25*L 0.805*L], 1, cs_wing; label = 2))
push!(meshes, frame_member([0 0.25*L 0.805*L;  0 8.5*L 0.805*L], nc, cs_wing; label = 2))
push!(meshes, frame_member([0 -0.25*L 0.805*L;  0 -8.5*L 0.805*L], nc, cs_wing; label = 2))
push!(meshes, frame_member([0 8.5*L 0.805*L;  0 9.5*L 0.805*L], 1, cs_wing; label = 2))
push!(meshes, frame_member([0 -8.5*L 0.805*L;  0 -9.5*L 0.805*L], 1, cs_wing; label = 2))
push!(meshes, frame_member([0 9.5*L 0.805*L;  0 10.0*L 0.805*L], 1, cs_wing; label = 2))
push!(meshes, frame_member([0 -9.5*L 0.805*L;  0 -10.0*L 0.805*L], 1, cs_wing; label = 2))

# Define the drums at the ends of the wings.
push!(meshes, frame_member([0 +9.5*L +0.91*L; +2*L +9.5*L +0.91*L], 1, cs_drum; label = 3))
push!(meshes, frame_member([0 +9.5*L +0.91*L; -2*L +9.5*L +0.91*L], 1, cs_drum; label = 3))
push!(meshes, frame_member([0 -9.5*L +0.91*L; +2*L -9.5*L +0.91*L], 1, cs_drum; label = 3))
push!(meshes, frame_member([0 -9.5*L +0.91*L; -2*L -9.5*L +0.91*L], 1, cs_drum; label = 3))

# Define the horizontal and vertical parts of the tail.
push!(meshes, frame_member([-8*L 0 .75*L; -8*L 0 3.35*L], 2, cs_tailv; label = 4))
push!(meshes, frame_member([-8*L 0 3.35*L; -8*L 0 3.75*L], 2, cs_tailv; label = 4))
push!(meshes, frame_member([-8*L 0 3.8*L; -8*L 2*L 3.8*L], 2, cs_tailh; label = 5))
push!(meshes, frame_member([-8*L 0 3.8*L; -8*L -2*L 3.8*L], 2, cs_tailh; label = 5))

# Define the parts of the aluminum constraining plate for the viscoelastic layer.
push!(meshes, frame_member([-.119*L 0 0.8665*L;  -.119*L 0.25*L 0.8665*L], 1, cs_vconstr; label = 6))
push!(meshes, frame_member([-.119*L 0 0.8665*L;  -.119*L -0.25*L 0.8665*L], 1, cs_vconstr; label = 6))
push!(meshes, frame_member([-.119*L 0.25*L 0.8665*L;  -.119*L 8.5*L 0.8665*L], nc, cs_vconstr; label = 6))
push!(meshes, frame_member([-.119*L -0.25*L 0.8665*L;  -.119*L -8.5*L 0.8665*L], nc, cs_vconstr; label = 6))


# Define the massless connectors between:
# Wing - Wingdrum
push!(meshes, frame_member([0 +9.5*L +0.805*L;0 +9.5*L +0.91*L], 1, cs_conn; label = 7))
push!(meshes, frame_member([0 -9.5*L +0.805*L;0 -9.5*L +0.91*L], 1, cs_conn; label = 7))
# Body-Wing
push!(meshes, frame_member([0 0 0; 0 0.0*L .805*L], 1, cs_connw; label = 7))
# Body-Tail
push!(meshes, frame_member([-8*L 0 0; -8*L 0 .75*L], 1, cs_connt; label = 7))
# Tail-Taildrum
push!(meshes, frame_member([-8*L 0 3.75*L; -8*L 0 3.8*L], 1, cs_conn; label = 7))

# Wing-Constraining plate for the viscoelastic layer
# Middle connector created individually
push!(meshes, frame_member([0 0 .805*L; -.119*L 0 0.8665*L], 1, cs_conn; label = 7))

# Connectors created for both wings
for i in 1:nc+1
    push!(meshes, frame_member([0 (0.25+(i-1)*8.25/nc)*L .805*L;  -.119*L (0.25+(i-1)*8.25/nc)*L 0.8665*L], 1, cs_conn; label = 7))
    push!(meshes, frame_member([0 -(0.25+(i-1)*8.25/nc)*L .805*L;  -.119*L -(0.25+(i-1)*8.25/nc)*L 0.8665*L], 1, cs_conn; label = 7))
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

# Wingdrum complementary masses
push!(meshes, frame_member([2*L 9.5*L .91*L ; 1.8*L 9.2*L .96*L], 1, cs_conn; label = 8))# added mass
push!(meshes, frame_member([2*L -9.5*L .91*L ; 1.8*L -9.2*L .96*L], 1, cs_conn; label = 8))# added mass

# Merge all the meshes of individual parts. This will glue together nodes which
# are in the "same" location. The parts of the mesh can be distinguished based
# on the label. 
fens, fesa = mergenmeshes(meshes, tolerance)

# The geometry is visualized in the tutorial [garteur_geometry_tut](garteur_geometry_tut.jl).

##
# ## Material

# Material properties can be now used to create a material: isotropic elasticity
# model of the `FinEtoolsDeforLinear` package is instantiated.
using FinEtoolsDeforLinear
# The material of the structure is aluminum.
# The elastic modulus:
E = 70000.0 * phun("MPa") 
nu = 0.31;
# The mass density:
rho = 2700 * phun("kg/m^3")
alu = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

# The material of the viscoelastic layer. The properties are due to the
# technical specification [3]. Since the properties are frequency dependent, we
# take as a representative value the numbers obtained for 20 Hz.
# Poisson ratio:
nu = 0.49;
# Storage modulus:
Gp = 0.6 * phun("MPa") 
E = 2 * (1 + nu) * Gp
# The mass density:
rho = 900 * phun("kg/m^3")
layer = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

# Material for the massless connectors has the mass density set to zero;
# otherwise it has the same properties as the aluminum material  of the
# structure.

massless = MatDeforElastIso(DeforModelRed3D, 0.0, alu.E, alu.nu, 0.0)

# This simple function returns material based on the label of the beam elements.
material(labl) = begin
    if labl == 6
        return layer
    elseif labl == 7 || labl == 8
        return massless
    end
    return alu
end

# This is the assumed stifffness of the bungee cords (each one separately).
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
# node. Note that the the incremental displacements are in general complex.
dchi = NodalField(0.0im .* zeros(size(fens.xyz, 1), 6))

# There are no support conditions.
applyebc!(dchi)
# The  the number of free
# (unknown) degrees of freedom is equal to the total number of degrees of freedom in the system.
numberdofs!(dchi);

##
# ## Identify support points and locations of sensors

# Suspension points
suspln = selectnode(fens; box = initbox!(Float64[], vec([0.0*L 0.0*L 0.805*L])), inflate = tolerance)
susprn = selectnode(fens; box = initbox!(Float64[], vec([0.0*L -0.0*L 0.805*L])), inflate = tolerance)
suspbn = selectnode(fens; box = initbox!(Float64[], vec([-2.0*L 0.0*L 0.0*L])), inflate = tolerance)


sensors = Dict()
sensors = let 
    # The sensors at the tip of the left and right wing drum
    sensor112n = selectnode(fens; box = initbox!(Float64[], vec([+1.8*L 9.8*L 0.96*L])), inflate = tolerance)
    sensors[112] = sensor112n
    sensor12n = selectnode(fens; box = initbox!(Float64[], vec([+1.8*L -9.8*L .96*L])), inflate = tolerance)
    sensors[12] = sensor12n
    sensor111n = selectnode(fens; box = initbox!(Float64[], vec([-1.8*L 9.8*L 0.96*L])), inflate = tolerance)
    sensors[111] = sensor111n
    sensor11n = selectnode(fens; box = initbox!(Float64[], vec([-1.8*L -9.8*L .96*L])), inflate = tolerance)
    sensors[11] = sensor11n
        # The joint between the horizontal and vertical tail parts
    sensor202n = selectnode(fens; box = initbox!(Float64[], vec([-8*L 0 3.8*L])), inflate = tolerance)
    sensors[202] = sensor202n
    sensors
end

##
# ## Assemble the global discrete system


# For disambiguation we will refer to the stiffness and mass functions by qualifying them with the corotational-beam module, `FEMMCorotBeamModule`.
using FinEtoolsFlexBeams.FEMMCorotBeamModule
CB = FEMMCorotBeamModule

# Note that we have an array of finite element sets. We compute the matrices for each set separately and accumulate them into the final overall matrix.
# Thus we can construct the stiffness and mass matrix as follows.
using  SparseArrays

Kf, Kd, M = let
    Kf = spzeros(dchi.nfreedofs, dchi.nfreedofs)
    Kd = spzeros(dchi.nfreedofs, dchi.nfreedofs)
    M = spzeros(dchi.nfreedofs, dchi.nfreedofs)
    for fes in fesa
        labl  = fes.label[1]
        femm = CB.FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material(labl));
        if labl == 6 # damping layer
            Kd += CB.stiffness(femm, geom0, u0, Rfield0, dchi);
        else
            Kf += CB.stiffness(femm, geom0, u0, Rfield0, dchi);
        end
        M += CB.mass(femm, geom0, u0, Rfield0, dchi);
    end
    Kf, Kd, M
end

##
# ## Additional concentrated masses.


using LinearAlgebra

using FinEtoolsFlexBeams.FEMMPointMassModule
PM = FEMMPointMassModule

# There is a sensor on the tail.
femmcm1 =  PM.FEMMPointMass(IntegDomain(FESetP1(reshape([sensor202n;], 1, 1)), PointRule()), FFltMat(2*L*L/5*L/5*2*rho*I(3)));

# These are the forward/interior locations on the wing drums.
mass1n = selectnode(fens; box = initbox!(Float64[], vec([1.8*L 9.2*L .96*L])), inflate = tolerance)
mass2n = selectnode(fens; box = initbox!(Float64[], vec([1.8*L -9.2*L .96*L])), inflate = tolerance)
femmcm2 =  PM.FEMMPointMass(IntegDomain(FESetP1(reshape([mass1n; mass2n;], 2, 1)), PointRule()), FFltMat(0.2*phun("kg")*I(3)));

Mp = PM.mass(femmcm1, geom0, u0, Rfield0, dchi) + PM.mass(femmcm2, geom0, u0, Rfield0, dchi);

##
# ## Bungee supports 


using LinearAlgebra

using FinEtoolsFlexBeams.FEMMPointGroundedSpringModule
BS = FEMMPointGroundedSpringModule

# There are three suspension points at the top of the fuselage. We assume that these bungee supports exert only reaction in the vertical direction.
femmbs =  BS.FEMMPointGroundedSpring(IntegDomain(FESetP1(reshape([suspln; susprn; suspbn;], 3, 1)), PointRule()), 
FFltMat([bungeecoefficient*[0;0;1]*[0;0;1]' 0*I(3); 0*I(3) 0*I(3)]));

Kb = BS.stiffness(femmbs, geom0, u0, Rfield0, dchi)



##
# ## Damping

# The measured overall damping ratios were approximately 1%. Damping levels
# thought appropriate for obtaining realistic damping were obtained through the
# use of a viscoelastic layer with an aluminum constraining layer. 



# ### Viscoelastic damping layer

# The viscoelastic used was the 3M acrylic viscoelastic polymer ISD 112 in the
# form of a 76 mm by 50 mm roll.  This material is particularly well suited for
# the testbed operating range of 5-50 Hz and at 20 degrees C where the loss
# factor is near its peak of 0.4 [3]. These are representative quantities taken
# at 20 Hz.

eta = 0.5
omega_5_50 = 2*pi*20
Cd = (eta / omega_5_50) .* Kd

# ### Aluminum structure

# We assume the loss factor of the fuselage, wing, and tail to be 0.01. This can
# be used to derive the damping model for these parts of the aircraft structure
# in the form of Rayleigh damping.
zeta= 0.01/2; # damping ratio
zeta_1= zeta; # damping ratio for mode 1
zeta_2= zeta; # damping ratio for mode 2
# Stiffness and mass proportional damping parameters.
omega_1=2*pi*6; # Guess
omega_2=2*pi*30; # Guess
rdmass = (2*1.0/(omega_2^2-omega_1^2)).*[omega_2^2 -omega_1^2]*[zeta_1*omega_1; zeta_2*omega_2];
rdstiffness = (2*1.0/(omega_2^2-omega_1^2)).*[-1 1]*[zeta_1*omega_1; zeta_2*omega_2];

Cf = rdmass .* M + rdstiffness .* Kf;

Kt = Kf + Kd + Kb
Mt = M + Mp
Ct = Cf + Cd

##
# ## Loading

# Here we assume that the stinger  was attached at the location of the sensor 12.
# The magnitude of the force is arbitrary.
forceat = 12
fmagn= 1.0;
loadbdry = FESetP1(reshape(sensors[forceat], 1, 1))
lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
# The force is applied in the vertical direction, and we assume it is positive
# upwards.
fi = ForceIntensity(FFlt[0, 0.0, -fmagn, 0, 0, 0]);

# The force loading is now integrated over the "volume" of the integration
# domain.
F = CB.distribloads(lfemm, geom0, dchi, fi, 3);

##
# ## Solve the Harmonic vibration problem

fromf = 3.0
tof = 50.0
nf = 150

frequencies = logspace(log10(fromf), log10(tof), nf);
receptance12 = fill(0.0im, nf)
receptance112 = fill(0.0im, nf)
v = fill(0.0im, 6)

for   fi in 1:length(frequencies)
    f =  frequencies[fi];
    om = 2*pi*f;
    U = (-om^2*Mt + 1im*om*Ct + Kt) \ F
    scattersysvec!(dchi, U)
    gathervalues_asvec!(dchi, v, sensors[12])
    p_d = v[1:3];
    receptance12[fi] = p_d[3]/fmagn;
    gathervalues_asvec!(dchi, v, sensors[112])
    p_d = v[1:3];
    receptance112[fi] = p_d[3]/fmagn;
end

oms = (2*pi) .* frequencies;
mobility12 = receptance12 .* (-1im*oms);
mobility112 = receptance112 .* (-1im*oms);
accelerance12 = receptance12 .* (-oms.^2);
accelerance112 = receptance112 .* (-oms.^2);

results = Dict()
results[12] = Dict("receptance"=>receptance12, "mobility"=>mobility12, "accelerance"=>accelerance12)
results[112] = Dict("receptance"=>receptance112, "mobility"=>mobility112, "accelerance"=>accelerance112)

# axes(controls.ax1);
# cla
# semilogy(frequencies,abs(H12),'bd-','linewidth',2,'markersize',2)
# hold on
# semilogy(frequencies,abs(H112),'rx-','linewidth',2,'markersize',2)
# xlabel('Frequency [Hz]');
# ylabel([type ' ' unit])
# grid on
# legend({'12','112'});

using PlotlyJS

quantity = "accelerance"; units = "m/s^2/N"
outputat = 12
y = abs.(results[outputat][quantity]) / phun(units)
tc12 = scatter(;x=frequencies, y=y, mode="lines", name = "output@$(outputat)", line_color = "rgb(15, 15, 15)")
outputat = 112
y = abs.(results[outputat][quantity]) / phun(units)
tc112 = scatter(;x=frequencies, y=y, mode="lines", name = "output@$(outputat)", line_color = "rgb(215, 15, 15)")
plots = cat(tc12, tc112; dims = 1)
layout = Layout(;width=650, height=400, xaxis=attr(title="Frequency [Hz]", type = "linear"), yaxis=attr(title="abs(H) [$(units)]", type = "log"), title = "Force@$(forceat), $(quantity)")
pl = plot(plots, layout; options = Dict(
        :showSendToCloud=>true, 
        :plotlyServerURL=>"https://chart-studio.plotly.com"
        ))
display(pl)

# axes(controls.ax2);
# cla
# plot(frequencies,atan2(imag(H12),real(H12))/pi*180,'bd-','linewidth',2,'markersize',2)
# hold on
# plot(frequencies,atan2(imag(H112),real(H112))/pi*180,'rx-','linewidth',2,'markersize',2)
# set(gca,'ylim',[-180,180])
# set(gca,'ytick',[-180:90:180])
# xlabel('Frequency [Hz]');
# ylabel('Phase [deg]')
# grid on
# legend({'12','112'});

true