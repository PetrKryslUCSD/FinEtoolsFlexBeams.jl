#
julia="/c/Users/PetrKrysl/AppData/Local/Programs/Julia 1.5.0-rc2/bin/julia"

for n in docs/src/tutorials/*.jl;                        
do          
    echo $(basename $n)     
    "$julia" -e "using Pkg; Pkg.activate(\".\"); Pkg.instantiate(); cd(\"docs/src/tutorials/\"); include(\"$(basename $n)\"); exit()"                
done      
