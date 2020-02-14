module RotUtilModule

using LinearAlgebra: norm
using FinEtools

function initial_Rfield(fens)
    f = NodalField(zeros(size(fens.xyz,1), 9))
    for i in 1:count(fens)
        f.values[i, 1] = f.values[i, 5] = f.values[i, 9] = 1.0
    end
    return f
end

end # module
