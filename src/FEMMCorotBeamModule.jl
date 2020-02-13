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
function _local_cartesian_to_natural(L)
    aN=[[ -1,   0,   0,  0,  0, 0, 1,    0,    0, 0,  0,  0]
        [  0,   0,   0,  0,  0, 1, 0,    0,    0, 0,  0, -1]
        [  0, 2/L,   0,  0,  0, 1, 0, -2/L,    0, 0,  0,  1]
        [  0,   0,   0,  0, -1, 0, 0,    0,    0, 0,  1,  0]
        [  0,   0, 2/L,  0, -1, 0, 0,    0, -2/L, 0, -1,  0]
        [  0,   0,   0, -1,  0, 0, 0,    0,    0, 1,  0,  0]];
end

# Compute forces through which the element acts on the nodes in the
# local coordinate system.
# 
function _local_forces(PN, L)
    aN= _local_cartesian_to_natural(L);
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
function _local_frames(x0,x1x2_vector,xt,RI,RJ)
    # This is the element frame in the configuration t=0
    F0 = fill(0.0, 3, 3);
    F0[:,1] = (x0(2,:)-x0[1,:])'; 
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
    Ft[:,1] = (xt[2,:]-xt[1,:])'; 
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
    LFtI=Ft'*FtI;
    LFtJ=Ft'*FtJ;
    
    # C  COMPUTE NET DEFORMATIONS- DELTA L,TH1IJ,TH2I,TH2J,TH3I,TH3J
    # C  CHANGE IN LENGTH
    dN = fill(0.0, 6);
    # Total change in length between configurations 0 and t
    dN[1] = Lt-L0;
    # Total axial TORSION ANGLE
    dN[6] = (LFtJ[3,2]/LFtJ[2,2]-LFtI[3,2]/LFtI[2,2] -LFtJ[2,3]/LFtJ[3,3]+LFtI[2,3]/LFtI[3,3])/2;
    # Total rotation  angles of nodal cross-sections relative to the element frame
    TH2I=-LFtI[3,1]/LFtI[1,1];
    TH2J=-LFtJ[3,1]/LFtJ[1,1];
    TH3I= LFtI[2,1]/LFtI[1,1];
    TH3J= LFtJ[2,1]/LFtJ[1,1];
    dN[2] = TH3I-TH3J; # symmetric bending
    dN[3] = TH3I+TH3J; # anti-symmetric bending
    dN[4] = -TH2I+TH2J; # symmetric bending
    dN[5] = -TH2I-TH2J; # anti-symmetric bending
    return Lt,Ft,dN
end


end # module
