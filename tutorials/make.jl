using Literate

for t in readdir(".")
    if occursin(r".*_tut.jl", t)
        println("\nTutorial $t in $(pwd())\n")
        Literate.markdown(t, "."; documenter=false);
        # Literate.notebook(t, "."; execute=false, documenter=false);
    end
end
