module curved_cantilever_examples
##  Curved (45 deg circular arc) cantilever, Transverse force at the tip
# 
# Large-deflection problem. Solved many times:
#
# Bathe, Boulourchi 1979
# Simo, Vu-Quoc 1986
# Cardona, Geradin 1988
# Krysl 1993
# Krysl, FAESOR script curved_cantilever
# 
# Present calculation refers to the data of
# Reference: EULERIAN FORMULATION FOR LARGE-DISPLACEMENT ANALYSIS OF SPACE
# FRAMES, by B.  A.  lzzuddin I and A.  S.  Elnashai,
# Journal of Engineering Mechanics, Vol.  119, No.  3,  March, 1993.


using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexBeams.FEMMCorotBeamModule
using FinEtoolsFlexBeams.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
geostiffness = FEMMCorotBeamModule.geostiffness
mass = FEMMCorotBeamModule.mass
distribloads_global = FEMMCorotBeamModule.distribloads_global
restoringforce = FEMMCorotBeamModule.restoringforce
qmass = FEMMCorotBeamModule.qmass
using FinEtoolsFlexBeams.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using FinEtoolsBeamsVis: plot_points, plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio
using PlotlyJS
using JSON


function curved_cantilever()
    # Parameters:
    E=1e7;# lb.in^-2
    nu=0.3;
    h= 1; b =1;# in
    Fmag=600;#  lb
    radius=100;# in
    ang=45;
    nel=8;# number of elements
    utol=1e-6;
    maxit = 150
    
    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0])
    
    # Select the number of elements per leg.
    fens, fes = frame_member([0 0 0; nel 0 0], nel, cs)
    for i in 1:count(fens)
        k = fens.xyz[i, 1]
        fens.xyz[i, :] .= (radius*cos(k/nel*ang*2*pi/360), radius*sin(k/nel*ang*2*pi/360), 0)
    end
    clampedn = [1]
    tipn = [nel+1];
    
    tolerance = radius/nel/100;
    
    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, 0.0, E, nu, 0.0)
    
    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))
    
    # Apply EBC's
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, clampedn, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    
    # Additional fields
    u1 = deepcopy(u0)
    Rfield1 = deepcopy(Rfield0)
    rhs = gathersysvec(dchi);
    TMPv = deepcopy(rhs)
    utol = 1e-13*dchi.nfreedofs;
    
    # tbox = plot_space_box([[0.0*radius 0.0*radius 0.0*radius]; [1.15*radius 0.75*radius 2.0*radius]])
    # tshape0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    # plots = cat(tbox,  tshape0; dims = 1)
    # pl = render(plots)
    # sleep(0.5)
    
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    loadbdry = FESetP1(reshape(tipn, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))
    
    load_parameters = 0.05:0.05:1;
    tip_displacement=[];
    for load_parameter in load_parameters
        applyebc!(dchi) # Apply boundary conditions
        u1.values[:] = u0.values[:]; # guess
        Rfield1.values[:] = Rfield0.values[:]; # guess
        fi = ForceIntensity(FFlt[0, 0, load_parameter*Fmag, 0, 0, 0]);
        
        println("Load: $load_parameter")
        iter = 1;
        while true
            F = distribloads(lfemm, geom0, dchi, fi, 3);
            Fr = restoringforce(femm, geom0, u1, Rfield1, dchi);       # Internal forces
            @. rhs = F + Fr;
            K = stiffness(femm, geom0, u1, Rfield1, dchi) + geostiffness(femm, geom0, u1, Rfield1, dchi);
            dchi = scattersysvec!(dchi, (K)\rhs); # Disp. incr
            u1.values[:] += (dchi.values[:,1:3])[:];   # increment displacement
            update_rotation_field!(Rfield1, dchi)
            @show iter, maximum(abs.(dchi.values[:]))
            if maximum(abs.(dchi.values[:])) < utol # convergence check
                break;
            end
            if (iter > maxit)# bailout for failed convergence
                error("Possible failed convergence");
            end
            iter += 1;
        end
        u0.values[:] = u1.values[:];       # update the displacement
        Rfield0.values[:] = Rfield1.values[:]; # update the rotations
        
        # if mod(step, 4) == 0
        #     @show step
        #     tshape1 = plot_solid(fens, fes; x = geom0.values, u = u1.values, R = Rfield1.values, facecolor = "rgb(125, 15, 15)");
        #     plots = cat(tbox,  tshape0,  tshape1; dims = 1)
        #     pl.plot.layout[:title] = "t=$t"
        #     react!(pl, plots, pl.plot.layout)
        #     sleep(0.5)
        # end
    end
    
    return true
end # simo_vuquoc_animated

function allrun()
    println("#####################################################")
    println("# slowtop1 ")
    slowtop1()
    return true
end # function allrun

end # module
