
using Pkg; Pkg.activate("."); Pkg.instantiate()


# include(".\\examples\\meshing\\mesh_examples.jl");   
# mesh_examples.curve_mesh_envelope()
# mesh_examples.argyr_l_frame_envelope()
# mesh_examples.argyr_l_frame_movable()

# include(".\\examples\\modal\\beam_modal_examples.jl"); 
# beam_modal_examples.allrun()


# include(".\\examples\\modal\\argyr_l_frame_modal_examples.jl");   
# argyr_l_frame_modal_examples.allrun()

include(".\\examples\\linear_buckling\\tippling_examples.jl");  
tippling_examples.tippling_1()