
using Pkg; Pkg.activate("."); Pkg.instantiate()


include(".\\examples\\meshing\\mesh_examples.jl");
mesh_examples.allrun()

include(".\\examples\\modal\\beam_modal_examples.jl");
beam_modal_examples.allrun()


include(".\\examples\\modal\\argyr_l_frame_modal_examples.jl");
# argyr_l_frame_modal_examples.allrun()
argyr_l_frame_modal_examples.argyr_l_frame_modal_anim()

include(".\\examples\\linear_buckling\\tippling_examples.jl");
tippling_examples.allrun()

include(".\\examples\\nonlin_transient\\fast_top_examples.jl");
fast_top_examples.fasttop1()
