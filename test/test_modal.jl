module mmodal1
using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexBeams.FEMMCorotBeamModule
using FinEtoolsFlexBeams.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
mass = FEMMCorotBeamModule.mass
using FinEtoolsFlexBeams.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using Test
function test()
        # Parameters:
        E=71240.0;#MPa
        nu=0.31;# Poisson ratio
        rho=5e-9;
        b=3.0; h=30.0; L=240.0; # cross-sectional dimensions and length of each leg in millimeters
        # Choose the mass formulation:
        mass_type=2;
        scale = 0.4

        # Reference frequencies
        reffs = [11.2732, 30.5269]
        neigvs = 2;

        # Cross-sectional properties
        cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 1.0, 0.0])

        # Select the number of elements per leg.
        xyz =
        n = 4;
        members = []
        push!(members, frame_member([0 0 L; L 0 L], n, cs))
        push!(members, frame_member([L 0 L; L 0 0], n, cs))
        fens, fes = merge_members(members; tolerance = L / 10000);

        # Material properties
        material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

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
        println("Natural frequencies: $fs [Hz]")
        println("Reference: $reffs [Hz]")

        @test norm(reffs - fs[1:length(reffs)]) ./ norm(reffs) <= 3.0e-5

        # Visualize vibration modes
        # mode = 2
        # tbox = plot_space_box([[-L/2 -L/2 0]; [L/2 L/2 1.1*L]])
        # tenv0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
        # plots = cat(tbox, tenv0; dims = 1)
        # pl = render(plots; title = "Mode $(mode)")
        # Rfield1 = deepcopy(Rfield0)
        # for xscale in scale.*sin.(collect(0:1:72).*(2*pi/21))
        #     scattersysvec!(dchi, xscale.*v[:, mode])
        #     Rfield1 = deepcopy(Rfield0)
        #     update_rotation_field!(Rfield1, dchi)
        #     tenv1 = plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield1.values, facecolor = "rgb(25, 255, 25)");
        #     plots = cat(tbox, tenv0, tenv1; dims = 1)
        #     react!(pl, plots, pl.plot.layout)
        #     sleep(0.115)
        # end

        return true
    end # argyr_l_frame_modal_anim

end
using .mmodal1
mmodal1.test()