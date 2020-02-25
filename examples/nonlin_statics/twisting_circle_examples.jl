module twisting_circle_examples
# Twisting of a circle into shape with two turns
# G. Rebel; Finite Rotation Shell Theory including Drill Rotations and its 
# Finite Element Implementation; PhD Thesis, Delft University Press (1998).

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
gyroscopic = FEMMCorotBeamModule.gyroscopic
using FinEtoolsFlexBeams.RotUtilModule: initial_Rfield, linear_update_rotation_field!, update_rotation_field!
using LinearAlgebra: dot
using Arpack
using LinearAlgebra
using SparseArrays
using FinEtoolsBeamsVis: plot_points, plot_nodes, plot_midline, render, plot_space_box, plot_solid, space_aspectratio
using PlotlyJS
using JSON
using ORCA

num3dig(a) = round(a * 1000) / 1000

function twisting_circle()
    # Parameters:
    E = 200000.0;
    nu = 0.3;# Poisson ratio
    rho = 7.8e-9;
    b = 0.6; h = 6.0;  # cross-sectional dimensions 
    radius = 120 # radius of the circle
    utol = 1e-3;
    maxit = 30
    rotmag = 2*pi
    # Select the number of elements around the circumference.
    nel = 64;
    
    # Cross-sectional properties
    cs = CrossSectionRectangle(s -> h, s -> b, s -> [0.0, 0.0, 1.0])

    tolerance = radius/nel/1000;
    fens, fes = frame_member([0 0 0; 2*pi 0 0], nel, cs)
    for i in 1:count(fens)
        a = fens.xyz[i, 1]
        fens.xyz[i, :] .= (radius+radius*cos(a), radius*sin(a), 0)
    end
    fens, fes = mergenodes(fens, fes, tolerance, [1, nel+1])
    
    # Material properties
    material = MatDeforElastIso(DeforModelRed3D, 0.0, E, nu, 0.0)
    
    # Construct the requisite fields, geometry and displacement
    # Initialize configuration variables
    geom0 = NodalField(fens.xyz)
    u0 = NodalField(zeros(size(fens.xyz,1), 3))
    Rfield0 = initial_Rfield(fens)
    dchi = NodalField(zeros(size(fens.xyz,1), 6))
    
    # Locate some nodes
    box = fill(0.0, 6)
    # The pinned end is
    initbox!(box, [0.0,0.0,0.0])
    clampedn = selectnode(fens; box = box, inflate = tolerance)
    # The force is applied at the elbow
    initbox!(box, [2*radius,0,0])
    torquen = selectnode(fens; box = box, inflate = tolerance)

    # Apply EBC's
    for i in [1, 2, 3, 4, 5, 6]
        setebc!(dchi, clampedn, true, i)
    end
    for i in [2, 3, 4, 5, 6]
        setebc!(dchi, torquen, true, i)
    end
    applyebc!(dchi)
    numberdofs!(dchi);
    
    # Additional fields
    u1 = deepcopy(u0)
    Rfield1 = deepcopy(Rfield0)
    rhs = gathersysvec(dchi);
    TMPv = deepcopy(rhs)
    utol = 1e-13*dchi.nfreedofs;
    
    femm = FEMMCorotBeam(IntegDomain(fes, GaussRule(1, 2)), material)
    
    tbox = plot_space_box([[0 -1.2*radius -1.2*radius]; [2.5*radius 1.2*radius 1.2*radius]])
    tshape0 = plot_solid(fens, fes; x = geom0.values, u = 0.0.*dchi.values[:, 1:3], R = Rfield0.values, facecolor = "rgb(125, 155, 125)", opacity = 0.3);
    plots = cat(tbox,  tshape0; dims = 1)
    pl = render(plots)
    sleep(0.5)
    savefig(pl, "twisting_circle-$(0).png", width = 300, height = 300)

    load_parameter_delta = 0.002
    load_parameters = 0.0:load_parameter_delta:1;
    step = 0
    for load_parameter in load_parameters
        dchi.values[:] .= 0.0
        dchi.fixed_values[torquen[1], 4] = rotmag * load_parameter_delta
        applyebc!(dchi) # Apply boundary conditions
        u1.values[:] = u0.values[:]; # guess
        Rfield1.values[:] = Rfield0.values[:]; # guess
        update_rotation_field!(Rfield1, dchi)

        println("Load: $load_parameter")
        iter = 1;
        while true
            Fr = restoringforce(femm, geom0, u1, Rfield1, dchi);       # Internal forces
            rhs = Fr;
            K = stiffness(femm, geom0, u1, Rfield1, dchi) + geostiffness(femm, geom0, u1, Rfield1, dchi);
            dchi = scattersysvec!(dchi, (K)\rhs); # Disp. incr
            dchi.values[torquen[1], 4] = 0.0
            u1.values[:] += (dchi.values[:,1:3])[:];   # increment displacement
            update_rotation_field!(Rfield1, dchi)
            print("$iter: ||du||=$(maximum(abs.(dchi.values[:])))\n")
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
        
        
        step = step + 1
        if mod(step, 1) == 0
            @show step
            tshape1 = plot_solid(fens, fes; x = geom0.values, u = u1.values, R = Rfield1.values, facecolor = "rgb(125, 15, 15)");
            plots = cat(tbox,  tshape0,  tshape1; dims = 1)
            pl.plot.layout[:title] = "Angle=$(num3dig(rotmag * load_parameter / pi))*pi"
            react!(pl, plots, pl.plot.layout)
            sleep(0.5)
            savefig(pl, "twisting_circle-$(step).png", width = 300, height = 300)
        end

    end
    step = step + 1
    
    plots = cat(
    scatter(;x=tip_displacement[:,1], y=load_parameters.*Fmag, mode="lines", line_color = "rgb(0, 0, 0)", name="Tip-x"),
    scatter(;x=tip_displacement[:,2], y=load_parameters.*Fmag, mode="lines", line_color = "rgb(0, 0, 0)", name="Tip-y"),
    scatter(;x=tip_displacement[:,3], y=load_parameters.*Fmag, mode="lines", line_color = "rgb(0, 0, 0)", name="Tip-z"),
    scatter(;x=izzelna_ux[:,1], y=izzelna_ux[:,2], mode="markers", line_color = "rgb(255, 0, 0)", name="Ref-x"),
    scatter(;x=izzelna_uy[:,1], y=izzelna_uy[:,2], mode="markers", line_color = "rgb(0, 255, 0)", name="Ref-y"),
    scatter(;x=izzelna_uz[:,1], y=izzelna_uz[:,2], mode="markers", line_color = "rgb(0, 0, 255)", name="Ref-z"); 
    dims = 1)
    layout = Layout(width=700, height=500, xaxis=attr(title="Tip displacement [in]"),
    yaxis=attr(title="Load [lb]"))
    pl = plot(plots, layout)
    display(pl)
    
    return true
end # curved_cantilever

function curved_cantilever_thin()
    # Parameters:
    E=1e7;# lb.in^-2
    nu=0.3;
    thinness_factor = 100
    h= 1.0/thinness_factor; b = 1.0/thinness_factor;# in
    Fmag=600.0/(thinness_factor^4);#  lb
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
    
    load_parameters = 0.0:0.05:1;
    tip_displacement = fill(0.0, length(load_parameters), 3);
    step = 0
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
            print("$iter: ||du||=$(maximum(abs.(dchi.values[:])))\n")
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
        
        
        step = step + 1
        tip_displacement[step, :] .= u1.values[tipn[1], :]

    end
    step = step + 1
    
    plots = cat(
    scatter(;x=tip_displacement[:,1], y=load_parameters.*Fmag, mode="lines", line_color = "rgb(0, 0, 0)", name="Tip-x"),
    scatter(;x=tip_displacement[:,2], y=load_parameters.*Fmag, mode="lines", line_color = "rgb(0, 0, 0)", name="Tip-y"),
    scatter(;x=tip_displacement[:,3], y=load_parameters.*Fmag, mode="lines", line_color = "rgb(0, 0, 0)", name="Tip-z"),
    scatter(;x=izzelna_ux[:,1], y=izzelna_ux[:,2]/(thinness_factor^4), mode="markers", line_color = "rgb(255, 0, 0)", name="Ref-x"),
    scatter(;x=izzelna_uy[:,1], y=izzelna_uy[:,2]/(thinness_factor^4), mode="markers", line_color = "rgb(0, 255, 0)", name="Ref-y"),
    scatter(;x=izzelna_uz[:,1], y=izzelna_uz[:,2]/(thinness_factor^4), mode="markers", line_color = "rgb(0, 0, 255)", name="Ref-z"); 
    dims = 1)
    layout = Layout(width=700, height=500, xaxis=attr(title="Tip displacement [in]"),
    yaxis=attr(title="Load [lb]"))
    pl = plot(plots, layout)
    display(pl)
    
    return true
end # curved_cantilever

function allrun()
    println("#####################################################")
    println("# slowtop1 ")
    slowtop1()
    return true
end # function allrun

end # module
