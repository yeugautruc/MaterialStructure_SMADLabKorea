clc;
L=0.5;
b=25e-3;
h=50e-3;
t=1e-3;
E=200e9;
V=0.3;
E1=E/(1-V^2);
G=76.9e9;
a=b*h*(h-b)/18;
c=48*t^3/(h+b);
b1=b*h*(b+h)*t/2;
b2=4*b^2*h^2*t/(h+b);
b4=b*h*(b-h)*t/2;
b5=2*b^2*h^2*t/(h+b);
b6=2*b^2*h^2*t/(h+b);
Ne= 5  ; 
%Number of elements
% element stiffness matrix
Kel=[   G*b1/L      -G*b4/2                 G*b5/L                     -G*b1/L       -G*b4/2                  -G*b5/L;
          -G*b4/2       G*b1*L/3+E*a/L    -G*b6/2                     G*b4/2       G*b1*L/6-E1*a/L      -G*b6/2;
           G*b5/L      -G*b6/2                 G*b2/L+E1*c*L/3      -G*b5/L      -G*b6/2                    -G*b2/L+E1*c*L/6;
           -G*b1/L      G*b4/2                 -G*b5/L                    -G*b1/L      G*b4/2                     G*b5/L;
           -G*b4/2      G*b1*L/6-E1*a/L   -G*b6/2                     G*b4/2      G*b1*L/3+E1*a/L       G*b6/2;
           -G*b5/L     -G*b6/2                 -G*b2/L+E1*c*L/6       G*b5/L       G*b6/2                     G*b2/L+E1*c*L/3 ];