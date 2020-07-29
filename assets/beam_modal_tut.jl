# # Modal analysis of beam

# Goals:
# - 

# Description:
# Vibration analysis of beam which is simply supported in one plane, and clamped in another. 
# The results are compared with analytical expressions.

# using LinearAlgebra: dot
# using Arpack
# using LinearAlgebra
# using SparseArrays
# using PlotlyJS

using FinEtools

# The material parameters may be defined with the specification of the units.
# The elastic properties are:
E = 30002.0*phun("ksi") 
nu = 0.0;
# The mass density is expressed using the gravitational acceleration and the weight of unit volume:
g = 32.17*12*phun("in/sec^2");
rho = 0.2*phun("lbm")/g;
# Here are the cross-sectional dimensions and the length of the beam between supports.
b = 18*phun("in"); h = 18*phun("in"); L = 300*phun("in");

# Analytical frequencies 
# 
# The analytical frequencies were taken from table 8-1 of Formulas for
# natural frequency and mode shape, Robert D. Blevins, Krieger publishing
# company, Malabar Florida, reprint edition 2001.
# 
# The beam is aligned with the ``Y`` global Cartesian coordinate. The beam
# behaves as a simply supported beam in the vertical plane (global Cartesian
# ``YZ``), while in the horizontal plane (global Cartesian ``XY``) it behaves
# as a clamped beam.
# 
# The cross-sectional properties are:
A = b*h;# in^2
I2 = b*h^3/12;# cm^4 
I3 = b^3*h/12;# cm^4 
# Then the analytical vibration frequencies for the first two modes are:
@show analyt_freq = [(1*pi)^2, (4.73004074)^2] .* (sqrt(E*I2/rho/A)/(2*pi*L^2));

# The purpose of the numerical model is to calculate approximation to these two
# analytical natural frequencies.
neigvs = length(analyt_freq);
 
# Cross-sectional properties are incorporated in the cross-section property. The
# three arguments supplied are functions. All are returning "constants". In
# particular the first two functions each return the dimension of the
# cross-section as a constant(the beam has a uniform cross-section); the third
# function defines the orientation of the cross-section in the global Cartesian
# coordinates. `[1.0, 0.0, 0.0]` is the vector that together with the tangent
# to the midline curve of the beam spans the ``x_1x_2`` plane of the local
# coordinates for the beam.
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionRectangle
cs = CrossSectionRectangle(s -> b, s -> h, s -> [1.0, 0.0, 0.0])

# Now we generate the mesh of the beam. The locations of its two endpoints are:
xyz = [[0 -L/2 0]; [0 L/2 0]]
# We will generate
n = 4
# beam elements along the member.
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member
fens, fes = frame_member(xyz, n, cs)
@show fens

using FinEtoolsFlexBeams.FEMMCorotBeamModule
using FinEtoolsFlexBeams.FEMMCorotBeamModule: FEMMCorotBeam
using FinEtoolsFlexBeams.FEMMCorotBeamModule: stiffness
using FinEtoolsFlexBeams.FEMMCorotBeamModule: mass
using FinEtoolsFlexBeams.RotUtilModule: initial_Rfield

# Material properties
material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)


##
# Construct the requisite fields, geometry and displacement
# Initialize configuration variables
geom0 = NodalField(fens.xyz)
u0 = NodalField(zeros(size(fens.xyz,1), 3))
Rfield0 = initial_Rfield(fens)
dchi = NodalField(zeros(size(fens.xyz,1), 6))

# Apply EBC's
l1 = selectnode(fens; box = [0 0 -L/2 -L/2 0 0], tolerance = L/10000)
for i in [1,2,3,5,6]
    setebc!(dchi, l1, true, i)
end
l1 = selectnode(fens; box = [0 0 L/2  L/2 0 0], tolerance = L/10000)
for i in [1,2,3,5,6]
    setebc!(dchi, l1, true, i)
end
applyebc!(dchi)
numberdofs!(dchi);

## 
# Choose the mass formulation:
mass_type=2;

# Assemble the global discrete system
femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
K = stiffness(femm, geom0, u0, Rfield0, dchi);
M = mass(femm, geom0, u0, Rfield0, dchi);

# Solve the eigenvalue problem
d,v,nev,nconv = eigs(K, M; nev=2*neigvs, which=:SM)
fs = real(sqrt.(complex(d)))/(2*pi)
println("Eigenvalues: $fs [Hz]")
  
# #
# for i=1:neigvs
#     disp(['  Eigenvector ' num2str(i) ' frequency ' num2str(Frequencies(i)) ' (analyt=' num2str(analyt_freq(i)) ') [Hz]']);
# end
# disp(' ')
return true
end # beam_modal

function allrun()
    println("#####################################################")
    println("# beam_modal ")
    beam_modal()
    return true
end # function allrun

end # module