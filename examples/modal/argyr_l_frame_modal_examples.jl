##  Modal analysis of the L-frame of Argyris.
# Eigenvector 1 frequency 11.2732  [Hz]
#   Eigenvector 2 frequency 30.5269  [Hz]
module argyr_l_frame_modal_examples

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexBeams.FEMMCorotBeamModule
using FinEtoolsFlexBeams.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
mass = FEMMCorotBeamModule.mass
using FinEtoolsFlexBeams.RotUtilModule: initial_Rfield
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using PlotlyJS
using Test

function argyr_l_frame_modal()
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
l1 = selectnode(fens; box = [0 0 0 0 L L], tolerance = L/10000)
for i in [1,2,3,4,5,6]
    setebc!(dchi, l1, true, i)
end
applyebc!(dchi)
numberdofs!(dchi);

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
#     disp(['  Eigenvector ' num2str(i) ' frequency ' num2str(Frequencies(i)) ' (analyt=' num2str(analyt_Frequencies(i)) ') [Hz]']);
# end
# disp(' ')
return true
end # argyr_l_frame_modal

function allrun()
    println("#####################################################")
    println("# argyr_l_frame_modal ")
    argyr_l_frame_modal()
    return true
end # function allrun

end # module