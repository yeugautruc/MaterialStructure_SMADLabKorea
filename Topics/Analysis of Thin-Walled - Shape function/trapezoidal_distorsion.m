clear;
clc;
syms a b h sina cosa real

k=[1 h/(2*sina) h^2/(4*sina^2) h^3/(8*sina^3) 0 0 0 0 0 0 0 0 0 0 0 0
    0 0 0 0 1 -a/2 a^2/4 -a^3/8 0 0 0 0 0 0 0 0
    0 0 0 0 1 a/2 a^2/4 a^3/8 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0 1 -h/(2*sina) h^2/(4*sina^2) -h^3/(8*sina^3) 0 0 0 0
    0 0 0 0 0 0 0 0 1 h/(2*sina) h^2/(4*sina^2) h^3/(8*sina^3) 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 0 0 1 -b/2 b^2/4 -b^3/8
    0 0 0 0 0 0 0 0 0 0 0 0 1 b/2 b^2/4 b^3/8
    1 -h/(2*sina) h^2/(4*sina^2) -h^3/(8*sina^3) 0 0 0 0 0 0 0 0 0 0 0 0
    0 1 h/(sina) (3/4)*h^2/sina^2 0 -1 a -(3/4)*a^2 0 0 0 0 0 0 0 0
    0 0 0 0 0 1 a (3/4)*a^2 0 -1 h/(sina) -(3/4)*h^2/sina^2 0 0 0 0
    0 0 0 0 0 0 0 0 0 1 h/(sina) (3/4)*h^2/sina^2  0 -1 b -(3/4)*b^2
    0 -1 h/(sina) -(3/4)*h^2/sina^2 0 0 0 0 0 0 0 0 0 1 b (3/4)*b^2
    0 0 2 3*h/sina 0 0 -2 3*a 0 0 0 0 0 0 0 0
    0 0 0 0 0 0 2 3*a 0 0 -2 3*h/sina 0 0 0 0
    0 0 0 0 0 0 0 0 0 0 2 3*h/sina 0 0 -2 3*b
    0 0 -2 3*h/sina 0 0 0 0 0 0 0 0 0 0 2 3*b];
f=[(cosa*h)/(a*sina)
    1
    1
    -(cosa*h)/(a*sina)
    -(cosa*h)/(b*sina)
    -1
    -1
    (cosa*h)/(b*sina)
    0
    0
    0
    0
    0
    0
    0
    0];
d=k\f