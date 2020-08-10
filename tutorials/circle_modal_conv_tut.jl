# # Modal analysis of a free-floating steel circle

# ## Description

# Vibration analysis of a beam simply supported in one plane, and clamped in
# another. The results are compared with analytical expressions. This is a
# benchmark from the NAFEMS Selected Benchmarks for Natural Frequency Analysis,
# publication: Test VM09: Circular Ring --  In-plane and Out-of-plane
# Vibration.

# ## Goals

# - Show convergence relative to reference values. 
# - Demonstrate the optimization of eigenvalue accuracy by choosing mass type.

# 

##
# ## Definition of the basic inputs

# The finite element code realize on the basic functionality implemented in this
# package.
using FinEtools

# The material parameters may be defined with the specification of the units.
# The elastic properties are:
E = 200.0 * phun("GPa") 
nu = 0.3;

# The mass density is
rho = 8000 * phun("kg/m^3")
# Here are the cross-sectional dimensions and the length of the beam between supports.
radius = 1.0 * phun("m"); diameter = 0.1 * phun("m"); 

##
# ## Reference frequencies 

# There will be 6 rigid body modes (zero natural frequencies).

# NAFEMS Finite Element Methods & Standards, Abbassian, F., Dawswell, D. J., and Knowles, N. C.
# Selected Benchmarks for Natural Frequency Analysis, Test No. 6. Glasgow: NAFEMS, Nov., 1987. 

# Mode                Reference Value (Hz)  NAFEMS Target Value (Hz)
# 7, 8 (out of plane)         51.85                 52.29 
# 9, 10 (in plane)            53.38                 53.97 
# 11, 12 (out of plane)      148.8                 149.7 
# 13, 14 (in plane)          151.0                 152.4 
# 15, 16 (out of plane)      287.0                 288.3 
# 17, 18 (in plane)          289.5                 288.3 

# The reference values are analytically determined (Blevins, FORMULAS FOR
# DYNAMICS, ACOUSTICS AND VIBRATION, Table 4.16). 

##
# ## Cross-section

# Cross-sectional properties are incorporated in the cross-section property. The
# three arguments supplied are functions. All are returning "constants". In
# particular the first two functions each return the dimension of the
# cross-section as a constant(the beam has a uniform cross-section); the third
# function defines the orientation of the cross-section in the global Cartesian
# coordinates. `[1.0, 0.0, 0.0]` is the vector that together with the tangent
# to the midline curve of the beam spans the $x_1x_2$ plane of the local
# coordinates for the beam.
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionCircle
cs = CrossSectionCircle(s -> diameter/2, s -> [1.0, 0.0, 0.0], 6/7) # Timoshenko
# cs = CrossSectionCircle(s -> diameter/2, s -> [1.0, 0.0, 0.0]) # Bernoulli
@show cs.parameters(0.0)

# The parameters of the structure:
R = radius
I = cs.parameters(0.0)[4]
m = rho * cs.parameters(0.0)[1]

# For instance the the first out of plane mode is listed in this table as
J = cs.parameters(0.0)[2]
G = E/2/(1+nu)
i = 2 # the first non-rigid body mode
@show i*(i^2-1)/(2*pi*R^2)*sqrt(E*I/m/(i^2+E*I/G/J))

# The first "ovaling" (in-plane) mode is:
i=2 # the first ovaling mode
@show i*(i^2-1)/(2*pi*R^2*(i^2+1)^(1/2))*sqrt(E*I/m)



# We will generate
# n = 20 # 52.2048, 52.2048, 53.7606, 53.7606, 148.8338
# n = 200 # 51.6087, 51.6087, 53.1098, 53.1098, 147.1189
# n = 2000 # 51.6028, 51.6028, 53.1033, 147.1022
n = 80
# beam elements along the member.
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member
tolerance = radius/n/1000;
fens, fes = frame_member([0 0 0; 2*pi 0 0], n, cs)
for i in 1:count(fens)
    a = fens.xyz[i, 1]
    fens.xyz[i, :] .= (radius+radius*cos(a), radius*sin(a), 0)
end
fens, fes = mergenodes(fens, fes, tolerance, [1, n+1])

##
# ## Material

using FinEtoolsDeforLinear
material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

##
# ## Fields

geom0 = NodalField(fens.xyz)
u0 = NodalField(zeros(size(fens.xyz, 1), 3))
using FinEtoolsFlexBeams.RotUtilModule: initial_Rfield
Rfield0 = initial_Rfield(fens)

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
# ## Assemble the global discrete system

using FinEtoolsFlexBeams.FEMMCorotBeamModule: FEMMCorotBeam
femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material);

# For disambiguation we will refer to the stiffness and mass functions by qualifying them with the corotational-beam module, `FEMMCorotBeamModule`.
using FinEtoolsFlexBeams.FEMMCorotBeamModule
CB = FEMMCorotBeamModule
K = CB.stiffness(femm, geom0, u0, Rfield0, dchi);
@show size(K)

##
# ## Solve the free-vibration problem

# The purpose of the numerical model is to calculate approximation to the reference frequencies.

neigvs = 18;

# The mass shift needs to be applied since the structure is free-floating.
oshift = (2*pi*15)^2

# The Arnoldi algorithm implemented in the well-known `Arpack` package is used
# to solve the generalized eigenvalue problem with the sparse matrices. As is
# common in structural dynamics, we request the smallest eigenvalues in
# absolute value (`:SM`). 
using Arpack

using FinEtoolsFlexBeams.FESetCorotBeamModule: MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA, 
MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA, 
MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA, 
MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA

results = let
    results = Dict()
    for mtype in [
        MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA, 
        MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA, 
        MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA, 
        MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA]
        M = CB.mass(femm, geom0, u0, Rfield0, dchi; mass_type = mtype);
        evals, evecs, nconv = eigs(K + oshift * M, M; nev=neigvs, which=:SM, ncv = 3*neigvs, maxiter = 2000);
        @assert nconv == neigvs
        results[mtype] = evals, evecs
    end
    results
end

sigdig(n) = round(n * 10000) / 10000

print("\nMASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA\n")
evals = results[MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA][1]
f = sqrt.([max(0, e - oshift) for e in evals]) / (2 * pi);
print("$(sigdig.(f))\n")

print("\nMASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA\n")
evals = results[MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA][1]
f = sqrt.([max(0, e - oshift) for e in evals]) / (2 * pi);
print("$(sigdig.(f))\n")

print("\nMASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA\n")
evals = results[MASS_TYPE_LUMPED_DIAGONAL_NO_ROTATION_INERTIA][1]
f = sqrt.([max(0, e - oshift) for e in evals]) / (2 * pi);
print("$(sigdig.(f))\n")

print("\nMASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA\n")
evals = results[MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA][1]
f = sqrt.([max(0, e - oshift) for e in evals]) / (2 * pi);
print("$(sigdig.(f))\n")

# 51.6425, 53.1497, 53.1497, 147.2318
# 51.6144, 53.1191, 53.1191, 147.1526
# 51.6074, 53.1115, 53.1115, 147.1328
using PlotlyJS
using FinEtools.AlgoBaseModule: richextrapol

# Modes 7 and 8
sols = [51.6425, 51.6144, 51.6074]
resextrap = richextrapol(sols, [4.0, 2.0, 1.0])  
errs = (sols .- resextrap[1])./resextrap[1]
t78 = scatter(;x=2*pi*radius./[80, 160, 320], y=errs, mode="markers+lines", name = "", line_color = "rgb(155, 15, 15)")

# Modes 9 and 10
sols = [53.1497, 53.1191, 53.1115]
resextrap = richextrapol(sols, [4.0, 2.0, 1.0])  
errs = (sols .- resextrap[1])./resextrap[1]
t910 = scatter(;x=2*pi*radius./[80, 160, 320], y=errs, mode="markers+lines", name = "", line_color = "rgb(15, 155, 15)")

# Modes 11 and 12
sols = [147.2318, 147.1526, 147.1328]
resextrap = richextrapol(sols, [4.0, 2.0, 1.0])  
errs = (sols .- resextrap[1])./resextrap[1]
t1112 = scatter(;x=2*pi*radius./[80, 160, 320], y=errs, mode="markers+lines", name = "", line_color = "rgb(15, 15, 155)")

layout = Layout(;width=650, height=400, xaxis=attr(title="Element size", type = "log"), yaxis=attr(title="Normalized error [ND]", type = "log"), title = "TB: Convergence of modes 7, ..., 12", xaxis_range=[-2, -1], yaxis_range=[-5, -2])
pl = plot([t78, t910, t1112], layout; options = Dict(
        :showSendToCloud=>true, 
        :plotlyServerURL=>"https://chart-studio.plotly.com"
        ))
display(pl)