module CrossSectionModule

using Dierckx
using LinearAlgebra: norm
using FinEtools

abstract  type AbstractCrossSectionType end

struct CrossSectionCircle{F} <: AbstractCrossSectionType
    shape::String
    # Function to return A, J, I1, I2, I3
    #     A = cross-sectional area of the rod
    #     I1= moment of inertia of the cross-section, twisting about the beam axis 
    #     I2, I3=moments of inertia, bending about the x2-axis (I2) or x3-axis (I3)
    #     J= torsion constant, (to give torsion stiffness =G*J)
    parameters::F
end

function CrossSectionCircle(radius)
    function parameters(s)
        R=radius(s);
        A=pi*R^2;
        J=pi/2*R^4;
        I1=pi/2*R^4;
        I2=pi/4*R^4;
        I3=pi/4*R^4;
        return A, J, I1, I2, I3
    end
    return CrossSectionCircle("circle", parameters)
end

struct CrossSectionHollowCircle{F} <: AbstractCrossSectionType
    shape::String
    parameters::F
end

function CrossSectionHollowCircle(innerradius, outerradius)
    function parameters(s)
        Rext=outerradius(s);
        Rint=innerradius(s);
        A=pi*(Rext^2-Rint^2);
        J=pi/2*(Rext^4-Rint^4);
        I2=pi/4*(Rext^4-Rint^4);
        I3=pi/4*(Rext^4-Rint^4);
        I1=I2+I3;
        J=pi/2*(Rext^4-Rint^4);
        return A, J, I1, I2, I3
    end
    return CrossSectionHollowCircle("hollow circle", parameters)
end

struct CrossSectionRectangle{F} <: AbstractCrossSectionType
    shape::String
    parameters::F
end

function CrossSectionRectangle(d2, d3)
    function parameters(s)
        d2=d2(s);
        d3=d3(s);
        a=max(d2,d3);
        b=min(d2,d3);
        itp=Spline1D([1,1.5,2.0, 2.5, 3.0, 4.0, 5.0, 6.0,10,20,40,80,200,2000], 
            [0.141, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.299, 0.312,0.317, 0.325, 0.33,1/3, 1/3]);
        c = itp(a/b)
        A=d2*d3;
        J=c*a*b^3;
        I2=d2*d3^3/12;
        I3=d2^3*d3/12;
        I1=I2+I3;
        return A, J, I1, I2, I3
    end
    return CrossSectionRectangle("rectangle", parameters)
end

end # module
