using Test
@time @testset "Section" begin
    include("test_section.jl")
end
@time @testset "Mesh" begin
    include("test_mesh.jl")
end
@time @testset "Modal" begin
    include("test_modal.jl")
end
@time @testset "Buckling" begin
    include("test_buckling.jl")
end
