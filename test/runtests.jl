module m1
using FinEtools
using FinEtoolsFrames: FESetL2CorotBeam
using Test
function test()
    fes = FESetL2CorotBeam([1 2; 2 3; 3 4])
    @show fes
end
end
using .m1
m1.test()

module m2
using FinEtools
using FinEtoolsFrames: FESetL2CorotBeam
using FinEtoolsFrames.MeshFrameMemberModule: frame_member
using GR
using Test
function test()
    L = 42
    xyz = [0 0 0;
    0 L/4 L*1/4;
    L/4 L/4 L*2/4;
    L/4 0 L*3/4;
    0 0 L]
   #  xyz = [      0         0         0
   # -2.1000    2.6880    2.1000
   # -2.9400    5.1240    4.2000
   # -2.7300    7.2660    6.3000
   # -1.6800    9.0720    8.4000
   #  0.0000   10.5000   10.5000
   #  2.1000   11.5080   12.6000
   #  4.4100   12.0540   14.7000
   #  6.7200   12.0960   16.8000
   #  8.8200   11.5920   18.9000
   # 10.5000   10.5000   21.0000
   # 11.5920    8.8200   23.1000
   # 12.0960    6.7200   25.2000
   # 12.0540    4.4100   27.3000
   # 11.5080    2.1000   29.4000
   # 10.5000         0   31.5000
   #  9.0720   -1.6800   33.6000
   #  7.2660   -2.7300   35.7000
   #  5.1240   -2.9400   37.8000
   #  2.6880   -2.1000   39.9000
   #       0         0   42.0000]
    nL = 20
    
    fens, fes = frame_member(xyz, nL; )
    @show fens, fes
    scatter3(fens.xyz[:, 1], fens.xyz[:, 2], fens.xyz[:, 3])
end
end
using .m2
m2.test()