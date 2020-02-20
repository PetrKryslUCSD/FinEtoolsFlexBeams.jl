module fast_top_examples
##  Fast Lagrangian top.
#
# Fast-spinning Lagrange  top simulated by a beam model. The reference
# solution is available in Lagrange_top_fast (rigid-body solver).

using FinEtools
using FinEtoolsDeforLinear
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionRectangle
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member, merge_members
using FinEtoolsFlexBeams.FEMMCorotBeamModule
using FinEtoolsFlexBeams.FEMMCorotBeamModule: FEMMCorotBeam
stiffness = FEMMCorotBeamModule.stiffness
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

function fasttop1()
    # Parameters:
    E = 71240.0;#MPa
    nu = 0.31;# Poisson ratio
    rho = 2.7e-9;
    Width = 60; Height = 60; Length = 4*Width;; # cross-sectional dimensions 
    Mass = Width*Height*Length*rho;
    IL =  (1/12)*Mass*(Width^2+Height^2);
    Ib =  (1/12)*Mass*(Length^2+Height^2)+Mass*(Length/2)^2;
    Ih =  (1/12)*Mass*(Width^2+Length^2)+Mass*(Length/2)^2;
    Omega0 = 313*pi;
    R0 = rotmat3([0.05 0 0]);
    g = 9.81e3; # 9.81 meters per second squared
    q = [0,0,-g*Width*Height*rho];
    utol = 1e-3;
    maxit = 12
    dt = min(2*pi/norm(Omega0)/10, 0.005);
    tend = 0.8;
    ng = 1/2; nb = 1/4*(1/2+ng)^2;
    # Choose the mass formulation:
    mass_type=2;

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> Width, s -> Width, s -> [1.0, 0.0, 0.0])

    # Select the number of elements per leg.
    spin_vector = R0*[0, 0, Omega0];
    X=[0 0 0;    reshape(R0*[0,0,Length], 1, 3)];
    n=8;
    tolerance=Length/n/100;
    members = []
    push!(members, frame_member(X, n, cs))
    fens, fes = merge_members(members; tolerance = tolerance);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

    
    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Apply EBC's
     box = fill(0.0, 6)
    initbox!(box, X[1,:])
    supportn = selectnode(fens; box = box, tolerance = tolerance)
    initbox!(box, X[2,:])
    tipn = selectnode(fens; box = box, tolerance = tolerance)
    for i in [1, 2, 3]
        setebc!(dchi, supportn, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);

    # Additional fields
    stepdchi = deepcopy(dchi)
    u1 = deepcopy(u0)
    v0 = deepcopy(dchi) # need the numbering of the degrees of freedom
    for i in 1:size(v0.values, 1)
        v0.values[i, 4:6] = spin_vector
    end
    v1 = deepcopy(dchi)
    a1 = deepcopy(dchi) # zero out the Acceleration
    a0 = deepcopy(dchi) # zero out the Acceleration
    Rfield1 = deepcopy(Rfield0)
    utol = 1e-13*dchi.nfreedofs;

    tbox = plot_space_box([[-1.1*Width -1.1*Width 0]; [1.1*Width 1.1*Width 1.1*Length]])
    tshape0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    plots = cat(tbox,  tshape0; dims = 1)
    pl = render(plots)
    sleep(3.5)

    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    fi = ForceIntensity(q);

    t = 0.0; #
    step = 0;
    while (t <= tend)
        t = t + dt;
        println("Time $(t)"); # pause
        # Initialization
        applyebc!(dchi) # Apply boundary conditions
        u1.values[:] = u0.values[:]; # guess
        Rfield1.values[:] = Rfield0.values[:]; # guess
        stepdchi.values[:] .= 0.0;# Total increment in current step
        a1.values[:] = -(1/nb/dt)*v0.values[:] -(1/2-nb)/nb*a0.values[:];
        v1.values[:] = v0.values[:] + dt*((1-ng)*a0.values[:] + ng*a1.values[:]);
        dchipv = dt*gathersysvec(v0) + (dt^2/2*(1-2*nb))*gathersysvec(a0)
        vpv = gathersysvec(v0) +(dt*(1-ng))*gathersysvec(a0);
        
        iter = 1;
        while true
            F = distribloads_global(femm, geom0, u1, Rfield1, dchi, fi)
            Fr = restoringforce(femm, geom0, u1, Rfield1, dchi);       # Internal forces
            F = F + Fr;
            K = stiffness(femm, geom0, u1, Rfield1, dchi);
            M = mass(femm, geom0, u1, Rfield1, dchi);
            G = qmass(femm, geom0, u1, Rfield1, v1, dchi);
            rhs = F + M*((-1/(nb*dt^2))*gathersysvec(stepdchi)+(1/(nb*dt^2))*dchipv) + G*((-ng/nb/dt)*gathersysvec(stepdchi)+(ng/nb/dt)*dchipv - vpv);
            dchi = scattersysvec!(dchi, (K+(ng/nb/dt)*G+(1/(nb*dt^2))*M)\rhs); # Disp. incr
            u1.values[:] += (dchi.values[:,1:3])[:];   # increment displacement
            stepdchi.values[:] += dchi.values[:]
            v1.values[:] += (ng/nb/dt)*dchi.values[:];
            a1.values[:] += (1/nb/dt^2)*dchi.values[:];
            update_rotation_field!(Rfield1, dchi)
            if maximum(abs.(dchi.values[:])) < utol# convergence check
                break; 
            end
            if (iter > maxit)# bailout for failed convergence
                error("Possible failed convergence");
            end
            iter += 1;
        end
        u0.values[:] = u1.values[:];       # update the displacement
        Rfield0.values[:] = Rfield1.values[:]; # update the rotations
        v0.values[:] = v1.values[:];       # update the velocities
        a0.values[:] = a1.values[:];       # update the accelerations
        
        tshape1 = plot_solid(fens, fes; x = geom0.values, u = u1.values, R = Rfield1.values, facecolor = "rgb(125, 15, 15)");
        plots = cat(tbox,  tshape0,  tshape1; dims = 1)
        react!(pl, plots, pl.plot.layout)
        sleep(0.115)

        step=step+1;
    end

    # # Visualize vibration modes
    # scattersysvec!(dchi, v[:, 1])
    # update_rotation_field!(Rfield0, dchi)
    # plots = cat(plot_nodes(fens),
    #     plot_solid(fens, fes; x = geom0.values, u = dchi.values[:, 1:3], R = Rfield0.values);
    #     dims = 1)
    # render(plots; aspectratio = space_aspectratio(fens.xyz))

    return true
end # fasttop1

function allrun()
    println("#####################################################")
    println("# fasttop1 ")
    fasttop1()
    return true
end # function allrun

end # module
