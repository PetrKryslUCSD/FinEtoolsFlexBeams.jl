# # Modal analysis of a toy airplane

# ## Description

# This virtual test application is based on the test article  
#     used by the GARTEUR Structures & Materials Action Group 19  
#     which organized a Round Robin exercise where 12 European laboratories  
#     tested a single structure between 1995 and 1997. The benchmark structure   
#     was a laboratory structure built to simulate the dynamic behaviour  
#     of an aeroplane. The structure was initially built for a benchmark  
#     study on experimental modal analysis conducted by the  
#     Structures and Materials Action Group (SM-AG19) of the Group  
#     for Aeronautical Research and Technology in EURope (GARTEUR).  
#         The test-bed was designed and manufactured by ONERA, France.

# ## Goals

# - Show how to construct model from multiple connected beams.

# 

##
# ## Definition of the basic inputs

# The finite element code realize on the basic functionality implemented in this
# package.
using FinEtools

# The material parameters may be defined with the specification of the units.
# The elastic properties are:
E = 70000.0 * phun("MPa") 
nu = 0.31;

# The mass density is
rho = 2700 * phun("kg/m^3")
# This is the characteristic length. The dimensions of the aircraft frame are
# expressed in terms of multiples of this characteristic unit.
L = 0.1*phun("m");


##
# ## Cross-section

# Cross-sectional properties are incorporated in the cross-section property. 
using FinEtoolsFlexBeams.CrossSectionModule: CrossSectionRectangle
# Body of the frame.
cs_body = CrossSectionRectangle(s -> 1.5*L, s -> L/2, s -> [1.0, 0.0, 1.0])
# Wing beam.
cs_wing = CrossSectionRectangle(s -> L/10, s -> L, s -> [0.0, 0.0, 1.0])
# Wing drums.
cs_drum = CrossSectionRectangle(s -> L/10, s -> L, s -> [0.0, 0.0, 1.0])
# Vertical part of the tail.
cs_tailv = CrossSectionRectangle(s -> L, s -> L/10, s -> [1.0, 0.0, 1.0])
# Horizontal part of the tail.
cs_tailh = CrossSectionRectangle(s -> L/10, s -> L, s -> [0.0, 0.0, 1.0])
# Massless connectors of the structural parts.
cs_conn = CrossSectionRectangle(s -> L/10, s -> L/10, s -> [1.0, 0.0, 1.0])
# Massless connector between the body and the wings.
cs_connw = CrossSectionRectangle(s -> L/2, s -> L/2, s -> [1.0, 0.0, 1.0])
# Massless connector between the body and the tail.
cs_connt = CrossSectionRectangle(s -> L, s -> L/5, s -> [1.0, 0.0, 1.0])
# Viscoelastic connecting layer.
cs_vconstr = CrossSectionRectangle(s -> L*(1.1/100), s -> L*(76.2/100), s -> [0.0, 0.0, 1.0])

# 
using FinEtoolsFlexBeams.MeshFrameMemberModule: frame_member
tolerance = L/10000;

# The parts of the mesh can be distinguished based on the label.
meshes = Tuple{FENodeSet, AbstractFESet}[]

# Define the constituent parts of the body of the aircraft.
push!(meshes, frame_member([-9*L 0 0; -8.5*L 0 0], 1, cs_body; label = 1))
push!(meshes, frame_member([-8.5*L 0 0; -8.0*L 0 0], 1, cs_body; label = 1))
push!(meshes, frame_member([-8.0*L 0 0; -2.0*L 0 0], 2, cs_body; label = 1))
push!(meshes, frame_member([-2.0*L 0 0; 0 0 0], 1, cs_body; label = 1))
push!(meshes, frame_member([0 0 0; 6*L 0 0], 2, cs_body; label = 1))

# Define the aluminum parts of the wings.
push!(meshes, frame_member([0 0 .805*L;  0 0.25*L .805*L], 1, cs_wing; label = 2))
push!(meshes, frame_member([0 0 .805*L;  0 -0.25*L .805*L], 1, cs_wing; label = 2))
for i in 1:15
    push!(meshes, frame_member([0 (0.25+(i-1)*5.75/15)*L .805*L;  0 (0.25+5.75/15*i)*L .805*L], 1, cs_wing; label = 2))
    push!(meshes, frame_member([0 -(0.25+(i-1)*5.75/15)*L .805*L;  0 -(0.25+5.75/15*i)*L .805*L], 1, cs_wing; label = 2))
end
for i in 1:8
    push!(meshes, frame_member([0 (6*L+(i-1)*.5*L) .805*L; 0 (6*L + 0.5*i*L) .805*L], 1, cs_wing; label = 2))
    push!(meshes, frame_member([0 (-6*L-(i-1)*.5*L) .805*L; 0 (-6*L - 0.5*i*L) .805*L], 1, cs_wing; label = 2))
end

# Define the drums at the ends of the wings.
push!(meshes, frame_member([0 +9.5*L +0.91*L; +2*L +9.5*L +0.91*L], 1, cs_drum; label = 3))
push!(meshes, frame_member([0 +9.5*L +0.91*L; -2*L +9.5*L +0.91*L], 1, cs_drum; label = 3))
push!(meshes, frame_member([0 -9.5*L +0.91*L; +2*L -9.5*L +0.91*L], 1, cs_drum; label = 3))
push!(meshes, frame_member([0 -9.5*L +0.91*L; -2*L -9.5*L +0.91*L], 1, cs_drum; label = 3))

# Define the horizontal and vertical parts of the tail.
push!(meshes, frame_member([-8*L 0 .75*L; -8*L 0 3.35*L], 2, cs_tailv; label = 4))
push!(meshes, frame_member([-8*L 0 3.35*L; -8*L 0 3.75*L], 2, cs_tailv; label = 4))
push!(meshes, frame_member([-8*L 0 3.8*L; -8*L 2*L 3.8*L], 2, cs_tailh; label = 5))
push!(meshes, frame_member([-8*L 0 3.8*L; -8*L -2*L 3.8*L], 2, cs_tailh; label = 5))

# Define the parts of the viscoelastic constraining layer.
push!(meshes, frame_member([-.119*L 0 .8606*L;  -.119*L 0.25*L .8606*L], 1, cs_vconstr; label = 6))
push!(meshes, frame_member([-.119*L 0 .8606*L;  -.119*L -0.25*L .8606*L], 1, cs_vconstr; label = 6))
for i in 1:15
    push!(meshes, frame_member([-.119*L (0.25+(i-1)*5.75/15)*L .8606*L;  -.119*L (0.25+5.75/15*i)*L .8606*L], 1, cs_vconstr; label = 6))
    push!(meshes, frame_member([-.119*L -(0.25+(i-1)*5.75/15)*L .8606*L;  -.119*L -(0.25+5.75/15*i)*L .8606*L], 1, cs_vconstr; label = 6))
end
for i in 1:5
    push!(meshes, frame_member([-.119*L (6*L+(i-1)*.5*L) .8606*L; -.119*L (6*L+0.5*i*L) .8606*L], 1, cs_vconstr; label = 6))
    push!(meshes, frame_member([-.119*L (-6*L-(i-1)*.5*L) .8606*L; -.119*L (-6*L-0.5*i*L) .8606*L], 1, cs_vconstr; label = 6))
end

# Define the massless connectors between:
# Wing - Wingdrum
push!(meshes, frame_member([0 +9.5*L +0.805*L;0 +9.5*L +0.91*L], 1, cs_conn; label = 7))
push!(meshes, frame_member([0 -9.5*L +0.805*L;0 -9.5*L +0.91*L], 1, cs_conn; label = 7))
# Body-Wing
push!(meshes, frame_member([0 0 0; 0 0 .805*L], 1, cs_connw; label = 7))
# Body-Tail
push!(meshes, frame_member([-8*L 0 0; -8*L 0 .75*L], 1, cs_connt; label = 7))
# Tail-Taildrum
push!(meshes, frame_member([-8*L 0 3.75*L; -8*L 0 3.8*L], 1, cs_conn; label = 7))

# Wing-Viscoelastic Constraining layer
# Middle connector created individually
push!(meshes, frame_member([0 0 .805*L; -.119*L 0 .8606*L], 1, cs_conn; label = 7))

# Connectors created for both wings
push!(meshes, frame_member([0 0.25*L .805*L;  -.119*L 0.25*L 0.8606*L], 1, cs_conn; label = 7))
push!(meshes, frame_member([0 0.25*L .805*L;  -.119*L -0.25*L 0.8606*L], 1, cs_conn; label = 7))
for i in 1:15
    push!(meshes, frame_member([0 (0.25+(i-1)*5.75/15)*L 0.805*L;  -.119*L (0.25+(i-1)*5.75/15)*L .8606*L], 1, cs_conn; label = 7))
    push!(meshes, frame_member([0 -(0.25+(i-1)*5.75/15)*L 0.805*L;  -.119*L -(0.25+(i-1)*5.75/15)*L .8606*L], 1, cs_conn; label = 7))
end
for i in 1:5
    push!(meshes, frame_member([0 (6*L+(i-1)*.5*L) .805*L; -.119*L (6*L+(i-1)*.5*L) .8606*L], 1, cs_conn; label = 7))
    push!(meshes, frame_member([0 (-6*L-(i-1)*.5*L) .805*L; -.119*L -(6*L+(i-1)*.5*L) .8606*L], 1, cs_conn; label = 7))
end

# Massless Sensor Connectors
# Tail Sensors
push!(meshes, frame_member([-8*L 2*L 3.8*L; -(153/20)*L (37/20)*L 3.85*L], 1, cs_conn; label = 8))# 303
push!(meshes, frame_member([-8*L -2*L 3.8*L; -(153/20)*L -(37/20)*L 3.85*L], 1, cs_conn; label = 8))# 301

# Wingdrum Sensors
push!(meshes, frame_member([0 9.5*L .91*L ; 0 9.8*L .96*L], 1, cs_conn; label = 8))# 101
push!(meshes, frame_member([-2*L 9.5*L .91*L ; -1.8*L 9.8*L .96*L], 1, cs_conn; label = 8))# 112
push!(meshes, frame_member([2*L 9.5*L .91*L ; 1.8*L 9.8*L .96*L], 1, cs_conn; label = 8))# 111

push!(meshes, frame_member([0 -9.5*L .91*L ; 0 -9.8*L .96*L], 1, cs_conn; label = 8))# 1
push!(meshes, frame_member([-2*L -9.5*L .91*L ; -1.8*L -9.8*L .96*L], 1, cs_conn; label = 8))# 12
push!(meshes, frame_member([2*L -9.5*L .91*L ; 1.8*L -9.8*L .96*L], 1, cs_conn; label = 8))# 11

# Merge all the meshes of individual parts. This will glued together notes which are in the "same" location.
fens, fesa = mergenmeshes(meshes, tolerance)

using PlotlyJS
using FinEtoolsFlexBeams.VisUtilModule: plot_solid, plot_space_box, render, default_layout_3d, save_to_json

colors = [
"rgb(125, 155, 155)",  # 1 body
"rgb(125, 155, 155)",  # 2 wing
"rgb(125, 155, 155)",  # 3 drums
"rgb(125, 155, 155)",  # 4 vertical tail
"rgb(125, 155, 155)",  # 5 horizontal tail
"rgb(125, 155, 15)",  # 6 viscoelastic constraining layer
"rgb(15, 15, 155)",  # 7 massless connectors
"rgb(125, 15, 15)",  # 8 sensor connectors
]

tbox = plot_space_box([[-1.2 * L -1.2 * L -1.2 * L]; [+1.2 * L +1.2 * L +1.2 * L]])
traces = let traces = tbox
    for fes in fesa
        labl  = fes.label[1]
        tm = plot_solid(fens, fes; facecolor=colors[labl]);
        traces = cat(traces, tm; dims = 1)
    end
    traces
end
layout = default_layout_3d(;width=900, height=900, options = Dict(:responsive=>true))
layout[:scene][:aspectmode] = "data"
pl = render(traces; layout = layout)

true