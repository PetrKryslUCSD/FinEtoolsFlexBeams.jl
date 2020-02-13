module FEMMCorotBeamModule

using LinearAlgebra: norm
using FinEtools

# Transformation from local Cartesian displacements to natural deformations 
# of a beam element. 
#
# Matrix defined in Eq (4.8) of COMPUTER  METHODS  IN APPLIED  MECHANICS  AND  ENGINEERING   14  (1978)  401-451
# ON LARGE  DISPLACEMENT-SMALL   STRAIN  ANALYSIS   OF  STRUCTURES
# WITH  ROTATIONAL   DEGREES  OF  FREEDOM.
# J.H.  ARGYRIS,   P.C.  DUNNE  and  D.W. SCHARPF
# 
# Inputs:
# L= current length of the element
# 
# Outputs:
# aN= transformation matrix to take Cartesian (local) displacement increments in the 
#      element frame and to produce increments of natural deformations; 
#      see local_frames() for the definition of the natural deformations 
function _local_cartesian_to_natural!(aN, L)
    fill!(aN, 0.0)
    aN[1, 1] = -1; aN[1, 7] = +1
    aN[2, 6] = +1; aN[2, 12] = -1
    aN[3, 2] = 2/L; aN[3, 6] = +1; aN[3, 8] = -2/L; aN[3, 12] = -1
    aN[4, 5] = -1; aN[4, 11] = +1
    aN[5, 3] = 2/L; aN[5, 5] = -1; aN[5, 9] = -2/L; aN[5, 11] = -1
    aN[6, 4] = -1; aN[6, 10] = +1
    # aN=[[ -1,   0,   0,  0,  0, 0, 1,    0,    0, 0,  0,  0]
    #     [  0,   0,   0,  0,  0, 1, 0,    0,    0, 0,  0, -1]
    #     [  0, 2/L,   0,  0,  0, 1, 0, -2/L,    0, 0,  0,  1]
    #     [  0,   0,   0,  0, -1, 0, 0,    0,    0, 0,  1,  0]
    #     [  0,   0, 2/L,  0, -1, 0, 0,    0, -2/L, 0, -1,  0]
    #     [  0,   0,   0, -1,  0, 0, 0,    0,    0, 1,  0,  0]];
        return aN
end

# Compute forces through which the element acts on the nodes in the
# local coordinate system.
# 
function _local_forces(PN, L, aN)
    _local_cartesian_to_natural(aN, L);
    FL=aN'*PN;
end

# Compute the current length of the element, the current element frame, 
# and the natural deformations  
#
# Inputs:
# x0= array of node coordinates, one node per row, in initial configuration
# x1x2_vector= vector that lies in the x1-x2 local element coordinate
#      plane, in initial configuration
# xt= array of node coordinates, one node per row, in initial configuration
# RI,RJ=nodal rotation (orthogonal) matrix 
# 
# Outputs:
# Lt= current length of the element,
# Ft= current element frame (orthogonal rotation matrix) whose columns are unit
#      vectors: they are centered halfway between the current locations of
#      the nodes, vector 1 points from node I to node J, vector 3 is
#      orthogonal to the sum of the nodal cross-section frame vectors 2
# dN= vector of natural deformations; dN(1)= total change in length 
#      between configurations 0 and t; dN(2)= symmetric bending;
#      dN(3)= anti-symmetric bending;    dN(4)= symmetric bending
#      dN(5)= anti-symmetric bending; dN(6)=total axial torsion angle.
# 
function _local_frames(x0, x1x2_vector, xt, RI, RJ)
    # This is the element frame in the configuration t=0
    F0 = fill(0.0, 3, 3);
    F0[:,1] = (x0(2,:)-x0[1,:]); 
    L0 = norm(@view F0[:,1]);
    F0[:,1] /= L0;
    #     F0(:,3)=skewmat(F0(:,1))*x1x2_vector;
    F0[:,3] .= (-F0[3,1]*x1x2_vector[2]+F0[2,1]*x1x2_vector[3],
                 F0[3,1]*x1x2_vector[1]-F0[1,1]*x1x2_vector[3],
                -F0[2,1]*x1x2_vector[1]+F0[1,1]*x1x2_vector[2]);
    F0[:,3] /= norm(@view F0[:,3]);
    #     F0(:,2)=skewmat(F0(:,3))*F0(:,1);
    F0[:,2] .= (-F0[3,3]*F0[2,1]+F0[2,3]*F0[3,1],
                 F0[3,3]*F0[1,1]-F0[1,3]*F0[3,1],
                -F0[2,3]*F0[1,1]+F0[1,3]*F0[2,1]);
             
    # The nodal cross-section frames are rotated by the nodal rotation matrices
    FtI = RI*F0;
    FtJ = RJ*F0;
    
    Ft = fill(0.0, 3, 3);
    # Compute the element frame in configuration t
    Ft[:,1] = (xt[2,:]-xt[1,:]); 
    Lt =norm(Ft(:,1));
    Ft[:,1] /= Lt;
    x1x2_vectort=FtI[:,2]+FtJ[:,2];
    #     Ft(:,3)=skewmat(Ft(:,1))*x1x2_vectort; Ft(:,3)=Ft(:,3)/norm(Ft(:,3));
    #     Ft(:,2)=skewmat(Ft(:,3))*Ft(:,1);
    #     Ft(:,3)=skewmat(Ft(:,1))*x1x2_vectort; # In the interest of speed,
    #     replace with below explicit rewrite
    Ft[:,3] = (-Ft[3,1]*x1x2_vectort[2]+Ft[2,1]*x1x2_vectort[3],
                Ft[3,1]*x1x2_vectort[1]-Ft[1,1]*x1x2_vectort[3],
               -Ft[2,1]*x1x2_vectort[1]+Ft[1,1]*x1x2_vectort[2];
    Ft[:,3] /= norm(Ft[:,3]);
    #     Ft(:,2)=skewmat(Ft(:,3))*Ft(:,1); # In the interest of speed,
    #     replace with below explicit rewrite
    Ft[:,2] = (-Ft[3,3]*Ft[2,1]+Ft[2,3]*Ft[3,1],
                Ft[3,3]*Ft[1,1]-Ft[1,3]*Ft[3,1],
               -Ft[2,3]*Ft[1,1]+Ft[1,3]*Ft[2,1]);
    
    # Components of FtI,FtJ in the element frame
    LFtI = Ft'*FtI;
    LFtJ = Ft'*FtJ;
    
    # C  COMPUTE NET DEFORMATIONS- DELTA L,TH1IJ,TH2I,TH2J,TH3I,TH3J
    # C  CHANGE IN LENGTH
    dN = fill(0.0, 6);
    # Total change in length between configurations 0 and t
    dN[1] = Lt-L0;
    # Total axial TORSION ANGLE
    dN[6] = (LFtJ[3,2]/LFtJ[2,2]-LFtI[3,2]/LFtI[2,2] -LFtJ[2,3]/LFtJ[3,3]+LFtI[2,3]/LFtI[3,3])/2;
    # Total rotation  angles of nodal cross-sections relative to the element frame
    TH2I =-LFtI[3,1]/LFtI[1,1];
    TH2J =-LFtJ[3,1]/LFtJ[1,1];
    TH3I = LFtI[2,1]/LFtI[1,1];
    TH3J = LFtJ[2,1]/LFtJ[1,1];
    dN[2] = TH3I-TH3J; # symmetric bending
    dN[3] = TH3I+TH3J; # anti-symmetric bending
    dN[4] = -TH2I+TH2J; # symmetric bending
    dN[5] = -TH2I-TH2J; # anti-symmetric bending
    return Lt,Ft,dN
end

# Compute the local geometric stiffness matrix. 
#
# function SM = local_geometric_stiffness(self,A,I2,I3,PN,L)
# 
# Inputs:
# A= cross-sectional area, 
# I2, I3=central moment of inertia of the cross-section about the x2 and x3 
# coordinate axis, 
# PN= vector of natural forces; see natural_forces() for definitions
# L= current length of the element, 
# 
# Outputs:
# SM = local geometric stiffness matrix, 12 x 12
# 
# This geometric stiffness matrix this consistent with relationship between
# the natural deformations and the natural forces that assumes there is
# only a linear constitutive link: no non-constitutive effects (bowing
# etc.) are included. This form of the geometric matrix was derived by
# Krenk.
# @BOOK{Krenk:2009,
#   AUTHOR =       {S. Krenk},
#   TITLE =        {Non-linear Modeling and Analysis of Solids and Structures },
#   PUBLISHER =    {Cambridge University Press},
#   YEAR =         {2009},
#   isbn =         {9780521830546}
# } 
function _local_geometric_stiffness!(SM, A, I2, I3, PN, L)
    N = PN[1];
    S_2 = -2*PN[3]/L;
    S_3 = -2*PN[5]/L;
    M_1 = PN[6];
    M_2I = PN[4]+PN[5];
    M_2J = PN[4]-PN[5];
    M_3I = -(PN[2]+PN[3]);
    M_3J = -(PN[2]-PN[3]);
    SM[1:3, 1:3] = [0   -S_2/L   -S_3/L;
                    -S_2/L   N/L   0;
                    -S_3/L   0   N/L];
    SM[1:3, 4:6] = [0   0   0;
                    -M_2I/L   M_1/L   0;
                    -M_3I/L   0   M_1/L]; 
    SM[1:3, 7:9] = [0   +S_2/L   +S_3/L;
                    +S_2/L   -N/L   0;
                    +S_3/L   0   -N/L]; 
    SM[1:3, 10:12] = [0   0   0;
                    M_2J/L   -M_1/L   0;
                    M_3J/L   0   -M_1/L]; 
    SM[4:6, 4:6] = [0   M_3I/2   -M_2I/2;
                    M_3I/2   0   0;
                    -M_2I/2   0   0];
    SM[4:6, 7:9] = [0   M_2I/L   M_3I/L;
                    0   -M_1/L   0;
                    0   0   -M_1/L]; 
    SM[4:6, 10:12] = [0   0   0;
                    0   0   M_1/2;
                    0   -M_1/2   0]; 
    SM[7:9, 7:9] = [0   -S_2/L   -S_3/L;
                    -S_2/L   N/L   0;
                    -S_3/L   0   N/L];
    SM[7:9, 10:12] = [0   0   0;
                    -M_2J/L   M_1/L   0;
                    -M_3J/L   0   M_1/L]; 
    SM[10:12, 10:12] = [0   -M_3J/2   M_2J/2;
                        -M_3J/2   0   0;
                        M_2J/2   0   0];
    complete_lt!(SM)
    return SM
end

# Mass matrix of the beam.
#
# function MM=local_mass (self, A, I1, I2, I3, rho, L)
#
# Inputs:
# A= cross-sectional area,
# I1=central moment of inertia of the cross-section about the x1 axis,
# I2, I3=central moment of inertia of the cross-section about the x2 and x3
# coordinate axis, 
# rho=mass density, 
# L0= initial length of the element, 
# 
# Outputs:
# MM = local mass matrix, 12 x 12
# In the element frame the mass matrix is constant.
function _local_mass!(MM, A, I1, I2, I3, rho, L, mass_type)
    if (mass_type == MASS_TYPE_CONSISTENT_WITH_ROTATION_INERTIA)
        # C
        # C  CONSISTENT MASS MATRIX including ROTATIONAL MASSES
        # C  Formulation of the (3.38), (3.39) equation from Dykstra's thesis
        MM .= (rho*A*L)*[...
            1/3      0        0        0        0        0       1/6       0        0        0        0        0
            0     13/35       0        0        0    11*L/210     0      9/70       0        0        0    -13*L/420
            0        0     13/35       0    -11*L/210    0        0        0       9/70      0    13*L/420     0
            0        0        0     I1/3/A      0        0        0        0        0      I1/6/A     0        0
            0        0   -11*L/210     0     L^2/105     0        0        0    -13*L/420    0    -L^2/140     0
            0    11*L/210     0        0        0     L^2/105     0    13*L/420     0        0        0    -L^2/140
            1/6      0        0        0        0        0       1/3       0        0        0        0        0
            0      9/70       0        0        0     13*L/420    0      13/35      0        0        0   -11*L/210
            0        0       9/70    0     -13*L/420     0        0        0      13/35      0    11*L/210     0
            0        0        0     I1/6/A      0        0        0        0        0      I1/3/A     0        0
            0        0    13*L/420     0    -L^2/140     0        0        0     11*L/210     0     L^2/105    0
            0   -13*L/420     0        0        0    -L^2/140     0   -11*L/210     0         0        0    L^2/105];
        MM .+= (rho/L)*[...
            0       0        0        0       0         0       0       0        0       0       0        0
            0    6/5*I2      0        0       0     L/10*I2     0    -6/5*I2     0       0       0    L/10*I2
            0       0     6/5*I3      0   -L/10*I3      0       0       0     -6/5*I3    0   -L/10*I3     0
            0       0        0        0       0         0       0       0        0       0       0        0
            0       0    -L/10*I3     0  2*L^2/15*I3    0       0       0     L/10*I3    0   -L^2/30*I3   0
            0    L/10*I2     0        0       0   2*L^2/15*I2   0    -L/10*I2    0       0       0    -L^2/30*I2
            0       0        0        0       0         0       0       0        0       0       0        0
            0   -6/5*I2      0        0       0     -L/10*I2    0    6/5*I2      0       0       0     -L/10*I2
            0       0     -6/5*I3     0    L/10*I3      0       0       0      6/5*I3    0    L/10*I3     0
            0       0        0        0       0         0       0       0        0       0       0        0
            0       0     -L/10*I3    0  -L^2/30*I3     0       0       0     L/10*I3    0   2*L^2/15*I3  0
            0   L/10*I2      0        0       0    -L^2/30*I2   0   -L/10*I2     0       0       0   2*L^2/15*I2];
    elseif (mass_type == MASS_TYPE_CONSISTENT_NO_ROTATION_INERTIA)
        # C
        # C  CONSISTENT MASS MATRIX excluding ROTATIONAL MASSES
        # C  Formulation of the (3.38) equation from Dykstra's thesis, no rotational inertia
        MM .= (rho*A*L)*[...
            1/3      0        0        0        0        0       1/6       0        0        0        0        0
            0     13/35       0        0        0    11*L/210     0      9/70       0        0        0    -13*L/420
            0        0     13/35       0    -11*L/210    0        0        0       9/70      0    13*L/420     0
            0        0        0     I1/3/A      0        0        0        0        0      I1/6/A     0        0
            0        0   -11*L/210     0     L^2/105     0        0        0    -13*L/420    0    -L^2/140     0
            0    11*L/210     0        0        0     L^2/105     0    13*L/420     0        0        0    -L^2/140
            1/6      0        0        0        0        0       1/3       0        0        0        0        0
            0      9/70       0        0        0     13*L/420    0      13/35      0        0        0   -11*L/210
            0        0       9/70    0     -13*L/420     0        0        0      13/35      0    11*L/210     0
            0        0        0     I1/6/A      0        0        0        0        0      I1/3/A     0        0
            0        0    13*L/420     0    -L^2/140     0        0        0     11*L/210     0     L^2/105    0
            0   -13*L/420     0        0        0    -L^2/140     0   -11*L/210     0         0        0    L^2/105];
    elseif (mass_type == MASS_TYPE_LUMPED_DIAGONAL_WITH_ROTATION_INERTIA)
        # C
        # C  LUMPED DIAGONAL MASS MATRIX WITH ROTATIONAL MASSES
        # C
        HLM  = A*rho*L/2.;
        HLI1 = rho*I1* L/2.;
        HLI2 = rho*I2* L/2.;
        HLI3 = rho*I3* L/2.;
        CA = HLM;
        CB = HLI1;
        CC = HLI2;
        CD = HLI3;
        fill!(MM, 0.0);
        MM[1,1]    = MM[1,1]    + CA;
        MM[2,2]    = MM[2,2]    + CA;
        MM[3,3]    = MM[3,3]    + CA;
        MM[4,4]    = MM[4,4]    + CB;
        MM[5,5]    = MM[5,5]    + CC;
        MM[6,6]    = MM[6,6]    + CD;
        MM[7,7]    = MM[7,7]    + CA;
        MM[8,8]    = MM[8,8]    + CA;
        MM[9,9]    = MM[9,9]    + CA;
        MM[10,10]  = MM[10,10]  + CB;
        MM[11,11]  = MM[11,11]  + CC;
        MM[12,12]  = MM[12,12]  + CD;
    elseif (self.mass_type == self.mass_type_lumped_diagonal_no_rotation_inertia)
        # C
        # C  LUMPED DIAGONAL ISOTROPIC MASS MATRIX WITHOUT ROTATIONAL MASSES
        # C
        HLM  = A*rho*L/2.;
        CA = HLM;
        CB = 0.0;
        CC = 0.0;
        CD = 0.0;
        fill!(MM, 0.0);
        MM[1,1]    = MM[1,1]    + CA;
        MM[2,2]    = MM[2,2]    + CA;
        MM[3,3]    = MM[3,3]    + CA;
        MM[7,7]    = MM[7,7]    + CA;
        MM[8,8]    = MM[8,8]    + CA;
        MM[9,9]    = MM[9,9]    + CA;
    end
    return MM
end


# Compute the local elastic stiffness matrix. 
#
# function SM = local_stiffness(self, E, G, A, I2, I3, J, L)
# 
# Inputs:
# E, G= Young's and shear modulus, 
# A= cross-sectional area, 
# I2, I3=central moment of inertia of the cross-section about the x2 and x3 
# coordinate axis, 
# J=St Venant torsion constant, 
# L= current length of the element, 
# 
# Outputs:
# SM = local stiffness matrix, 12 x 12
function _local_stiffness(SM, E, G, A, I2, I3, J, L, aN)
    _local_cartesian_to_natural!(aN, L);
    _natural_stiffness!(DN, E, G, A, I2, I3, J, L);
    SM .= aN'*DN*aN;
    return SM
end

# Compute the natural forces from the natural deformations.
#
# function PN = natural_forces (self, E, G, A, I2, I3, J, L, dN)
#
# Inputs:
# E, G= Young's and shear modulus, 
# A= cross-sectional area, 
# I2, I3=central moment of inertia of the cross-section about the x2 and x3 
# coordinate axis, 
# J=St Venant torsion constant, 
# L= current length of the element, 
# dN= column vector of natural deformations; see local_frames() 
# 
# Outputs:
#    PN = column vector of natural forces; 
#      PN(1)= axial force; 
#      PN(2)= symmetric bending moment in the plane x1-x2; 
#      PN(3)= anti-symmetric bending bending moment in the plane x1-x2;    
#      PN(4)= symmetric bending bending moment in the plane x1-x3;
#      PN(5)= anti-symmetric bending bending moment in the plane x1-x3; 
#      PN(6)= axial torque. 
# 
function _natural_forces!(PN, E, G, A, I2, I3, J, L, dN, DN)
    _natural_stiffness!(DN, E, G, A, I2, I3, J, L);
    #     Natural forces
    PN .= DN*dN;
    # Note that the non-constitutive stiffness due to pre-existing internal forces is currently omitted
end

# Compute the natural stiffness matrix.
#
# function DN= natural_stiffness(self, E, G, A, I2, I3, J, L) 
#
# Inputs:
# E, G= Young's and shear modulus, 
# A= cross-sectional area, 
# I2, I3=central moment of inertia of the cross-section about the x2 and x3 
# coordinate axis, 
# J=St Venant torsion constant, 
# L= current length of the element, 
# 
# Outputs:
# DN = 6 x 6 natural stiffness matrix
# 
function _natural_stiffness!(DN, E, G, A, I2, I3, J, L)
    fill!(DN, 0.0)
    DN[1, 1] = E*A
    DN[2, 2] = E*I3,
    DN[3, 3] = 3*E*I3,
    DN[4, 4] = E*I2,
    DN[5, 5] = 3*E*I2,
    DN[6, 6] = G*J ]/L;
    return DN
end

end # module
