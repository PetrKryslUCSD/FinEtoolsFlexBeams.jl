module CrossSectionModule

using Dierckx
using LinearAlgebra: norm
using FinEtools

abstract  type AbstractCrossSectionType end

struct CrossSectionCircle{F} <: AbstractCrossSectionType
    shape::String
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

# % Calculate the cross section parameters for a given shape of the cross section.
# %
# % function P=Beam_cross_section(shape) 
# % 
# % Input: 
# % shape = structure with attributes
# %      type= Type of the cross-section (as string), either
# %           circle: Additional attributes R (radius of the circle)
# %           hollow_circle: Additional attributes Rint,Rext 
# %                     (Internal and external radius of the circle) 
# %           square: Additional attributes a (length of the side)
# %           rectangle: Additional attributes d2,d3 (lengths of the sides)
# % 
# % Example:
# % P=Beam_cross_section(struct( 'type','hollow_circle','Rint',r,'Rext',r+t) ) 
# % P=Beam_cross_section(struct( 'type','rectangle','d2',b,'d3',h) ) 
# % 
# % Output:
# % P = structure with attributes 
# %     type = the same string as above
# %     A = cross-sectional area of the rod
# %     I1= moment of inertia of the cross-section, twisting about the beam axis 
# %     I2, I3=moments of inertia, bending about the x2-axis (I2) or x3-axis (I3)
# %     J= torsion constant, (to give torsion stiffness =G*J)
# function P=Beam_cross_section(shape)
#     P.type=shape.type;
#     if (isfield(shape,'x1x2_vector'))
#         P.x1x2_vector=shape.x1x2_vector;
#     end
#     switch shape.type
#         case 'circle'
#             P.R=shape.R;
#             P.A=pi*shape.R^2;
#             P.J=pi/2*shape.R^4;
#             P.I1=pi/2*shape.R^4;
#             P.I2=pi/4*shape.R^4;
#             P.I3=pi/4*shape.R^4;
#             P.external_radius=shape.R;
#         case 'hollow_circle'
#             P.Rext=shape.Rext;
#             P.Rint=shape.Rint;
#             P.A=pi*(shape.Rext^2-shape.Rint^2);
#             P.J=pi/2*(shape.Rext^4-shape.Rint^4);
#             P.I2=pi/4*(shape.Rext^4-shape.Rint^4);
#             P.I3=pi/4*(shape.Rext^4-shape.Rint^4);
#             P.I1=P.I2+P.I3;
#             P.J=pi/2*(shape.Rext^4-shape.Rint^4);
#             P.external_radius=shape.Rext;
#         case 'square'
#             P.a=shape.a;
#             P.A=shape.a^2;
#             P.J=0.1406*shape.a^4;
#             P.I1=shape.a^4/6;
#             P.I2=shape.a^4/12;
#             P.I3=shape.a^4/12;
#             P.external_radius=sqrt(2)/2*shape.a;
#         case 'rectangle'
#             P.d2=shape.d2;
#             P.d3=shape.d3;
#             a=max([shape.d2,shape.d3]);
#             b=min([shape.d2,shape.d3]);
#             c=interp1([1,1.5,2.0, 2.5, 3.0, 4.0, 5.0, 6.0,10,20,40,80,200,2000],...
#                 [0.141, 0.196, 0.229, 0.249, 0.263, 0.281, 0.291, 0.299, 0.312,0.317, 0.325, 0.33,1/3, 1/3],a/b);
#             P.A=shape.d2*shape.d3;
#             P.J=c*a*b^3;
#             P.I2=shape.d2*shape.d3^3/12;
#             P.I3=shape.d2^3*shape.d3/12;
#             P.I1=P.I2+P.I3;
#             P.external_radius=sqrt(shape.d2^2+shape.d3^2)/2;
#         otherwise
#             error([' Unknown type: ' char(shape.type)]);
#     end
# end


end # module
