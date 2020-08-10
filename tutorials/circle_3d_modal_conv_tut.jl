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


# The parameters of the structure:
# R = radius
# I = cs.parameters(0.0)[4]
# m = rho * cs.parameters(0.0)[1]

# # For instance the the first out of plane mode is listed in this table as
# J = cs.parameters(0.0)[2]
# G = E/2/(1+nu)
# i = 2 # the first non-rigid body mode
# @show i*(i^2-1)/(2*pi*R^2)*sqrt(E*I/m/(i^2+E*I/G/J))

# # The first "ovaling" (in-plane) mode is:
# i=2 # the first ovaling mode
# @show i*(i^2-1)/(2*pi*R^2*(i^2+1)^(1/2))*sqrt(E*I/m)



# We will generate
nperradius, nL = 8, 320
# nperradius, nL = 4, 160
# nperradius, nL = 2, 80
tolerance = diameter/nperradius/100
# beam elements along the member.
fens, fes = H8cylindern(diameter/2, 2*pi, nperradius, nL)
z0l = selectnode(fens, plane = [0, 0, 1, 0], inflate = tolerance)
z2l = selectnode(fens, plane = [0, 0, 1, 2*pi], inflate = tolerance)

for i in 1:count(fens)
    a = fens.xyz[i, 3]
    x = fens.xyz[i, 1]
    y = fens.xyz[i, 2] + radius
    fens.xyz[i, :] .= (x, y*cos(a), y*sin(a))
end
fens, fes = mergenodes(fens, fes, tolerance, vcat(z0l, z2l))


# File = "ring.vtk"
# vtkexportmesh(File, fens, fes)
# @async run(`"paraview.exe" $File`)

##
# ## Material

using FinEtoolsDeforLinear
MR = DeforModelRed3D

material = MatDeforElastIso(MR, rho, E, nu, 0.0)

femm = FEMMDeforLinearESNICEH8(MR, IntegDomain(fes, NodalTensorProductRule(3)), material)

geom = NodalField(fens.xyz)
u = NodalField(zeros(size(fens.xyz,1),3)) # displacement field

using SymRCM
p = let
    C = connectionmatrix(femm, count(fens))
    p = symrcm(C)
end

numberdofs!(u, p)

associategeometry!(femm,  geom)
K  = stiffness(femm, geom, u)
M = mass(femm, geom, u)

neigvs = 18;

# The mass shift needs to be applied since the structure is free-floating.
oshift = (2*pi*15)^2

using Arpack

evals, evecs, nconv = eigs(K+oshift*M, M; nev=neigvs, which=:SM)
evals = evals .- oshift;
sigdig(n) = round(n * 10000) / 10000
fs = real(sqrt.(complex(evals)))/(2*pi)
println("Eigenvalues: $(sigdig.(fs)) [Hz]")

# 2, 80, 59.7122, 63.142, 170.771
# 4, 160 53.7307, 55.8224, 153.3047
# 8, 320, 52.0987, 53.879

using FinEtools.AlgoBaseModule: richextrapol
richextrapol([59.7122,  53.7307, 52.0987], [4.0, 2.0, 1.0])  
richextrapol([63.142, 55.8224, 53.879], [4.0, 2.0, 1.0])  
richextrapol([170.771, 153.3047, 148.5574], [4.0, 2.0, 1.0])  

using PlotlyJS
using FinEtools.AlgoBaseModule: richextrapol

# Modes 7 and 8
sols = [59.7122,  53.7307, 52.0987]
resextrap = richextrapol(sols, [4.0, 2.0, 1.0])  
errs = (sols .- resextrap[1])./resextrap[1]
t78 = scatter(;x=2*pi*radius./[80, 160, 320], y=errs, mode="markers+lines", name = "Mode 7,8", line_color = "rgb(155, 15, 15)")

# Modes 9 and 10
sols = [63.142, 55.8224, 53.879]
resextrap = richextrapol(sols, [4.0, 2.0, 1.0])  
errs = (sols .- resextrap[1])./resextrap[1]
t910 = scatter(;x=2*pi*radius./[80, 160, 320], y=errs, mode="markers+lines", name = "Mode 9,10", line_color = "rgb(15, 155, 15)")

# Modes 11 and 12
sols = [170.771, 153.3047, 148.5574]
resextrap = richextrapol(sols, [4.0, 2.0, 1.0])  
errs = (sols .- resextrap[1])./resextrap[1]
t1112 = scatter(;x=2*pi*radius./[80, 160, 320], y=errs, mode="markers+lines", name = "Mode 11, 12", line_color = "rgb(15, 15, 155)")

layout = Layout(;width=400, height=300, xaxis=attr(title="Element size", type = "log"), yaxis=attr(title="Normalized error [ND]", type = "log"), title = "3D: Convergence of modes 7, ..., 12", xaxis_range=[-2, -1], yaxis_range=[-3, -0])
pl = plot([t78, t910, t1112], layout; options = Dict(
        :showSendToCloud=>true, 
        :plotlyServerURL=>"https://chart-studio.plotly.com"
        ))
display(pl)