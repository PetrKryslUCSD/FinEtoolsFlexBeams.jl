module argyr_swing_examples
# Argyris L-frame, large-amplitude swinging motion.
#
# An L-shaped frame is pushed at the elbow by a triangular-pulse force,
# which sends it swinging about the pinned end. The frame executes
# large-amplitude dynamic motion.

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

function argyr_swing_1()
    # Parameters:
    E = 71240.0;#MPa
    nu = 0.31;# Poisson ratio
    rho = 5.0e-9;
    b=10; h=30; L=240; #cross-sectional dimensions and length of each leg in millimeters
    utol = 1e-3;
    maxit = 30
    Fmag = 500;
    dt = 0.001;
    tend = 0.5;
    ng = 1/2; nb = 1/4*(1/2+ng)^2;
    # Choose the mass formulation:
    mass_type=1;

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0])

    # Select the number of elements per leg.
    n=1;
    tolerance=L/n/100;
    members = []
    push!(members, frame_member([0 L 0;L L 0], n, cs))
    push!(members, frame_member([L L 0;L 0 0], n, cs))
    fens, fes = merge_members(members; tolerance = tolerance);
    fens.xyz = Float64[240   240     0
            240     0     0
            0   240     0]
    fes.conn = [(3, 1), (1, 2)]
    @show fes
    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Locate some nodes
    box = fill(0.0, 6)
    # The pinned end is
    initbox!(box, [0,L,0])
    @show pinnedn = selectnode(fens; box = box, tolerance = tolerance)
    # The force is applied at the elbow
    initbox!(box, [L,L,0])
    @show elbown = selectnode(fens; box = box, tolerance = tolerance)
    # Tip is tracked
    initbox!(box, [L,0,0])
    @show tipn = selectnode(fens; box = box, tolerance = tolerance)


    # Apply EBC's
    for i in [1, 2, 3]
        setebc!(dchi, pinnedn, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    dchi.dofnums = [     1     2     3     4     5     6
     7     8     9    10    11    12
     0     0     0    13    14    15]

    # Additional fields
    stepdchi = deepcopy(dchi)
    u1 = deepcopy(u0)
    v0 = deepcopy(dchi) # need the numbering of the degrees of freedom
    v1 = deepcopy(dchi)
    a1 = deepcopy(dchi) # zero out the Acceleration
    a0 = deepcopy(dchi) # zero out the Acceleration
    Rfield1 = deepcopy(Rfield0)
    v0v = gathersysvec(v0)
    a0v = gathersysvec(a0)
    vpv = gathersysvec(v0);
    dchipv = gathersysvec(dchi);
    stepdchiv = gathersysvec(dchi);
    rhs = gathersysvec(dchi);
    TMPv = deepcopy(rhs)
    utol = 1e-13*dchi.nfreedofs;

    tbox = plot_space_box([[-1.1*L -1.1*L -1.1*L]; [1.1*L 1.1*L 1.1*L]])
    tshape0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    plots = cat(tbox,  tshape0; dims = 1)
    pl = render(plots)
    sleep(0.5)

    # tipx = Float64[]
    # tipy = Float64[]
    # push!(tipx, X[2,1])
    # push!(tipy, X[2,2])
    # tbox = scatter(;x=[-1.0, 1.0], y=[-1.0, 1.0], mode="markers")
    # plots = cat(tbox, scatter(;x=tipx./Length, y=tipy./Length, mode="markers+lines"); dims = 1)
    # layout = Layout(width=600, height=500)
    # pl = plot(plots, layout)
    # display(pl)

    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    loadbdry = FESetP1(reshape(elbown, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))

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
        gathersysvec!(v0, v0v)
        gathersysvec!(a0, a0v)
        dchipv = dt*v0v + (dt^2/2*(1-2*nb))*a0v
        vpv = v0v +(dt*(1-ng))*a0v;
        fi = ForceIntensity(FFlt[0, 0, (t<0.2)*((t<0.1)*t+(t>0.1)*(0.2-t))*Fmag, 0, 0, 0]);

        iter = 1;
        while true
            F = distribloads(lfemm, geom0, dchi, fi, 3);
            Fr = restoringforce(femm, geom0, u1, Rfield1, dchi);       # Internal forces
            @. rhs = F + Fr;
            K = stiffness(femm, geom0, u1, Rfield1, dchi);
            M = mass(femm, geom0, u1, Rfield1, dchi);
            G = qmass(femm, geom0, u1, Rfield1, v1, dchi);
            gathersysvec!(stepdchi, stepdchiv)
            @. TMPv = ((-1/(nb*dt^2))*stepdchiv+(1/(nb*dt^2))*dchipv)
            rhs .+= M*TMPv
            @. TMPv = ((-ng/nb/dt)*stepdchiv+(ng/nb/dt)*dchipv - vpv)
            rhs .+= G*TMPv;
             dchi = scattersysvec!(dchi, (K+(ng/nb/dt)*G+(1/(nb*dt^2))*M)\rhs); # Disp. incr
            u1.values[:] += (dchi.values[:,1:3])[:];   # increment displacement
            stepdchi.values[:] += dchi.values[:]
            v1.values[:] += (ng/nb/dt)*dchi.values[:];
            a1.values[:] += (1/nb/dt^2)*dchi.values[:];
            update_rotation_field!(Rfield1, dchi)
            @show maximum(abs.(dchi.values[:]))
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
        v0.values[:] = v1.values[:];       # update the velocities
        a0.values[:] = a1.values[:];       # update the accelerations

        if (mod(step,2)==0)
            tshape1 = plot_solid(fens, fes; x = geom0.values, u = u1.values, R = Rfield1.values, facecolor = "rgb(125, 15, 15)");
            plots = cat(tbox,  tshape0,  tshape1; dims = 1)
            react!(pl, plots, pl.plot.layout)
        end

        # if (mod(step,20)==0)
        #     push!(tipx, X[2,1]+u1.values[tipn[1], 1])
        #     push!(tipy, X[2,2]+u1.values[tipn[1], 2])
        #     plots = cat(tbox, scatter(;x=tipx./Length, y=tipy./Length, mode="markers+lines"); dims = 1)
        #     react!(pl, plots, pl.plot.layout)
        #     sleep(0.01)
        # end

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
end # argyr_swing_animated

function argyr_swing_animated()
    # Parameters:
    E = 71240.0;#MPa
    nu = 0.31;# Poisson ratio
    rho = 5.0e-9;
    b=10; h=30; L=240; #cross-sectional dimensions and length of each leg in millimeters
    utol = 1e-3;
    maxit = 30
    Fmag = 500;
    dt = 0.001;
    tend = 0.5;
    ng = 1/2; nb = 1/4*(1/2+ng)^2;
    # Choose the mass formulation:
    mass_type=1;

    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> b, s -> h, s -> [0.0, 0.0, 1.0])

    # Select the number of elements per leg.
    n=4;
    tolerance=L/n/100;
    members = []
    push!(members, frame_member([0 L 0;L L 0], n, cs))
    push!(members, frame_member([L L 0;L 0 0], n, cs))
    fens, fes = merge_members(members; tolerance = tolerance);

    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, rho, E, nu, 0.0)

    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))

    # Locate some nodes
    box = fill(0.0, 6)
    # The pinned end is
    initbox!(box, [0,L,0])
    pinnedn = selectnode(fens; box = box, tolerance = tolerance)
    # The force is applied at the elbow
    initbox!(box, [L,L,0])
    elbown = selectnode(fens; box = box, tolerance = tolerance)
    # Tip is tracked
    initbox!(box, [L,0,0])
    tipn = selectnode(fens; box = box, tolerance = tolerance)


    # Apply EBC's
    for i in [1, 2, 3]
        setebc!(dchi, pinnedn, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    dchi.dofnums

    # Additional fields
    stepdchi = deepcopy(dchi)
    u1 = deepcopy(u0)
    v0 = deepcopy(dchi) # need the numbering of the degrees of freedom
    v1 = deepcopy(dchi)
    a1 = deepcopy(dchi) # zero out the Acceleration
    a0 = deepcopy(dchi) # zero out the Acceleration
    Rfield1 = deepcopy(Rfield0)
    v0v = gathersysvec(v0)
    a0v = gathersysvec(a0)
    vpv = gathersysvec(v0);
    dchipv = gathersysvec(dchi);
    stepdchiv = gathersysvec(dchi);
    rhs = gathersysvec(dchi);
    TMPv = deepcopy(rhs)
    utol = 1e-13*dchi.nfreedofs;

    tbox = plot_space_box([[-1.1*L -0.1*L -1.1*L]; [1.1*L 2.1*L 1.1*L]])
    tshape0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    plots = cat(tbox,  tshape0; dims = 1)
    pl = render(plots)
    sleep(0.5)

    # tipx = Float64[]
    # tipy = Float64[]
    # push!(tipx, X[2,1])
    # push!(tipy, X[2,2])
    # tbox = scatter(;x=[-1.0, 1.0], y=[-1.0, 1.0], mode="markers")
    # plots = cat(tbox, scatter(;x=tipx./Length, y=tipy./Length, mode="markers+lines"); dims = 1)
    # layout = Layout(width=600, height=500)
    # pl = plot(plots, layout)
    # display(pl)

    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    loadbdry = FESetP1(reshape(elbown, 1, 1))
    lfemm = FEMMBase(IntegDomain(loadbdry, PointRule()))

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
        gathersysvec!(v0, v0v)
        gathersysvec!(a0, a0v)
        dchipv = dt*v0v + (dt^2/2*(1-2*nb))*a0v
        vpv = v0v +(dt*(1-ng))*a0v;
        fi = ForceIntensity(FFlt[0, 0, (t<0.2)*((t<0.1)*t+(t>0.1)*(0.2-t))*Fmag, 0, 0, 0]);

        iter = 1;
        while true
            F = distribloads(lfemm, geom0, dchi, fi, 3);
            Fr = restoringforce(femm, geom0, u1, Rfield1, dchi);       # Internal forces
            @. rhs = F + Fr;
            K = stiffness(femm, geom0, u1, Rfield1, dchi);
            M = mass(femm, geom0, u1, Rfield1, dchi);
            G = qmass(femm, geom0, u1, Rfield1, v1, dchi);
            gathersysvec!(stepdchi, stepdchiv)
            @. TMPv = ((-1/(nb*dt^2))*stepdchiv+(1/(nb*dt^2))*dchipv)
            rhs .+= M*TMPv
            @. TMPv = ((-ng/nb/dt)*stepdchiv+(ng/nb/dt)*dchipv - vpv)
            rhs .+= G*TMPv;
            dchi = scattersysvec!(dchi, (K+(ng/nb/dt)*G+(1/(nb*dt^2))*M)\rhs); # Disp. incr
            u1.values[:] += (dchi.values[:,1:3])[:];   # increment displacement
            stepdchi.values[:] += dchi.values[:]
            v1.values[:] += (ng/nb/dt)*dchi.values[:];
            a1.values[:] += (1/nb/dt^2)*dchi.values[:];
            update_rotation_field!(Rfield1, dchi)
            @show maximum(abs.(dchi.values[:]))
            if true
                tshape1 = plot_solid(fens, fes; x = geom0.values, u = u1.values, R = Rfield1.values, facecolor = "rgb(125, 155, 15)");
                plots = cat(tbox,  tshape0,  tshape1; dims = 1)
                react!(pl, plots, pl.plot.layout)
            end
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
        v0.values[:] = v1.values[:];       # update the velocities
        a0.values[:] = a1.values[:];       # update the accelerations

        if (mod(step,2)==0)
            tshape1 = plot_solid(fens, fes; x = geom0.values, u = u1.values, R = Rfield1.values, facecolor = "rgb(125, 15, 15)");
            plots = cat(tbox,  tshape0,  tshape1; dims = 1)
            react!(pl, plots, pl.plot.layout)
            sleep(0.1)
        end

        # if (mod(step,20)==0)
        #     push!(tipx, X[2,1]+u1.values[tipn[1], 1])
        #     push!(tipy, X[2,2]+u1.values[tipn[1], 2])
        #     plots = cat(tbox, scatter(;x=tipx./Length, y=tipy./Length, mode="markers+lines"); dims = 1)
        #     react!(pl, plots, pl.plot.layout)
        #     sleep(0.01)
        # end

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
end # argyr_swing_animated

function allrun()
    println("#####################################################")
    println("# slowtop1 ")
    slowtop1()
    return true
end # function allrun

end # module
