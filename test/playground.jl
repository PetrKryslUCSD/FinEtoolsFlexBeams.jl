module mmassform13
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexBeams.FEMMCorotBeamModule
using FinEtoolsFlexBeams.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
mass = FEMMCorotBeamModule.mass
geostiffness = FEMMCorotBeamModule.geostiffness
using FinEtoolsFlexBeams.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using InteractiveUtils
using Traceur
using FinEtoolsFlexBeams.FESetCorotBeamModule: FESetL2CorotBeam, local_frame_and_def!, local_mass!, local_stiffness!, natural_forces!, local_geometric_stiffness!, local_forces!
using FinEtoolsFlexBeams.FESetCorotBeamModule: local_mass_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA!, MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA

using Test
function test()
    # Parameters:
    E = 1000000.0;
    nu = 0.3;
    rho = 2700
    L1 = L0 = L =   30.0; # Length of the beam
    b =  0.5; # width
    h =  4.0; # height
    ecoords0 = fill(0.0, 2, 3); 
    ecoords1 = fill(0.0, 2, 3)
    edisp1 = fill(0.0, 2, 3); 
    evel1 = fill(0.0, 2, 6); 
    evel1f = fill(0.0, 2, 6)
    dofnums = zeros(FInt, 1, 12); 
    F0 = fill(0.0, 3, 3); 
    Ft = fill(0.0, 3, 3); 
    FtI = fill(0.0, 3, 3); 
    FtJ = fill(0.0, 3, 3)
    Te = fill(0.0, 12, 12)
    tempelmat1 = fill(0.0, 12, 12); 
    tempelmat2 = fill(0.0, 12, 12); 
    tempelmat3 = fill(0.0, 12, 12)
    elmat = fill(0.0, 12, 12);    
    elmatTe = fill(0.0, 12, 12);    
    elmato = fill(0.0, 12, 12)
    elvec = fill(0.0, 12);    
    elvecf = fill(0.0, 12)
    aN = fill(0.0, 6, 12)
    dN = fill(0.0, 6)
    DN = fill(0.0, 6, 6)
    PN = fill(0.0, 6)
    LF = fill(0.0, 12)
    R1I = fill(0.0, 3, 3);    
    R1J = fill(0.0, 3, 3);    
    OS = fill(0.0, 3, 3)

    # Reference frequencies
    reffs = [48.5475, 124.1839]

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [-1.0, 0.0, 0.0])

    # Select the number of elements per leg.
    n = 32;
    members = []
    push!(members, frame_member([0 0 0; 0 0 L], n, cs))
    fens, fes = merge_members(members; tolerance = L / 10000);

    A, I1, I2, I3, J, x1x2_vector = fes.A, fes.I1, fes.I2, fes.I3, fes.J, fes.x1x2_vector
    G = E / 2 / (1 + nu)
    i = 1

    A, I1, I2, I3, rho, L = 1.3, 35.1, 32.3, 16.5, 10000.0, 3.333
    MM = fill(0.0, 12, 12)
    local_mass_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA!(MM, A, I1, I2, I3, rho, L)
    MMp = fill(0.0, 12, 12)
    local_mass!(MMp, A, I1, I2, I3, rho, L, MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA)
    @test norm(MM - MMp) / norm(MM) <= 1.0e-6
    true
end
end
using .mmassform13
mmassform13.test()
