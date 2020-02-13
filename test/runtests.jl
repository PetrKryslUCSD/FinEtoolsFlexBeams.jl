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

module mcrosssection1
using FinEtools
using FinEtoolsFrames.CrossSectionModule
using Test
function test()
    cs = CrossSectionModule.CrossSectionCircle(s -> 5.9910)
    @show A, J, I1, I2, I3 = cs.parameters(0.0)
    for (c, r) in zip((A, J, I1, I2, I3), (112.75829799164978, 2023.5649824692157, 2023.5649824692157, 1011.7824912346078, 1011.7824912346078))
       @test c ≈ r
   end 
   true
end
end
using .mcrosssection1
mcrosssection1.test()

module mcrosssection2
using FinEtools
using FinEtoolsFrames.CrossSectionModule
using Test
function test()
    cs = CrossSectionModule.CrossSectionRectangle(s -> 42.0, s -> 4.2)
    @show A, J, I1, I2, I3 = cs.parameters(0.0)
    for (c, r) in zip((A, J, I1, I2, I3), (176.4, 970.849152, 26190.108000000004, 259.30800000000005, 25930.800000000003))
       @test c ≈ r
   end 
   true
end
end
using .mcrosssection2
mcrosssection2.test()

module mcrosssection3
using FinEtools
using FinEtoolsFrames.CrossSectionModule
using Test
function test()
    cs = CrossSectionModule.CrossSectionRectangle(s -> 1.3 * 4.2, s -> 4.2)
    @show A, J, I1, I2, I3 = cs.parameters(0.0)
    for (c, r) in zip((A, J, I1, I2, I3), (22.932000000000006, 71.66370760423138, 90.68000760000004, 33.71004000000001, 56.969967600000025))
       @test c ≈ r
   end 
   true
end
end
using .mcrosssection3
mcrosssection3.test()