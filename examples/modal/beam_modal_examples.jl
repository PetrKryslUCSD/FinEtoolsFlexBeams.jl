##  Simply-supported/clamped beam modal analysis. Corotational beam
# 
# Vibration analysis of beam which is simply supported in one plane, and clamped in another. 
# The results are compared with analytical expressions.
module beam_modal_examples

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFrames.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFrames.MeshFrameMemberModule: frame_member
using FinEtoolsFrames.FEMMCorotBeamModule
using FinEtoolsFrames.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using PlotlyJS
using Test

function beam_modal()
## 
# Parameters:
E=30002.0*1000;#30002ksi
nu=0.0;
b=18; h=18; L=300;# in, cross-section and length
g=32.17*12; # in/second^2
rho=0.2/g;# Pound/g mass density per unit volume

## 
# Cross-sectional properties
cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0])

## 
# Choose the mass formulation:
mass_type=2;

## 
# Analytical frequencies 
# 
# The analytical frequencies were taken from table 8-1 of Formulas for
# natural frequency and mode shape, Robert D. Blevins, Krieger publishing
# company, Malabar Florida, reprint edition 2001.
# 
# The beam has cylindrical supports at either end. Therefore in the
# vertical plane it behaves as a simply supported beam, while in the
# horizontal plane it behaves as a clamped beam.
# 
# Vertical plane, simply supported:
# Horizontal plane, clamped:
A=b*h;# in^2
I2=b*h^3/12;# cm^4 
I3=b^3*h/12;# cm^4 
analyt_Frequencies = [(1*pi)^2/(2*pi*L^2)*sqrt(E*I2/rho/A),
(4.73004074)^2/(2*pi*L^2)*sqrt(E*I2/rho/A)];
neigvs = length(analyt_Frequencies);
 
# Select the number of elements per half the length.
xyz = [[0 -L/2 0];[0 L/2 0]]
n=4;
fens, fes = frame_member(xyz, n, cs)

## 
# Material properties
material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)

##
# Construct the requisite fields, geometry and displacement
# Initialize configuration variables
geom0 = NodalField(fens.xyz)
u0 = NodalField(zeros(size(fens.xyz,1), 3))
Rfield0 = NodalField(zeros(size(fens.xyz,1), 9))
dchi = NodalField(zeros(size(fens.xyz,1), 6))

# Apply EBC's
l1 = selectnode(fens; box = [0 -L/2 0 0 -L/2 0], tolerance = L/10000)
for i in [1,2,3,5,6]
    setebc!(dchi, l1, true, i)
end
l1 = selectnode(fens; box = [0 L/2 0 0 L/2 0], tolerance = L/10000)
for i in [1,2,3,5,6]
    setebc!(dchi, l1, true, i)
end
applyebc!(dchi)
numberdofs!(dchi);

# Assemble the global discrete system
Km = stiffness(femm, geom0, u0, Rfield0, dchi);
# Mm = mass(femm,sysmat_assembler_sparse,geom0,u0,Rfield0,dchi);

# Solve the eigenvalue problem
# d,v,nev,nconv = eigs(K+OmegaShift*M, M; nev=neigvs, which=:SM)
# d = d .- OmegaShift;
# fs = real(sqrt.(complex(d)))/(2*pi)
# println("Eigenvalues: $fs [Hz]")
  
# #
# for i=1:neigvs
#     disp(['  Eigenvector ' num2str(i) ' frequency ' num2str(Frequencies(i)) ' (analyt=' num2str(analyt_Frequencies(i)) ') [Hz]']);
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