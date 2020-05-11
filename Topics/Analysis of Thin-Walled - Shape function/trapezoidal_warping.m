clear;
clc;
syms h sina a b cosa real

k=[h/(2*sina) 1 a/2 -1 0 0 0 0
    0 0 a/2 1 h/(2*sina) -1 0 0
    0 0 0 0 h/(2*sina) 1 b/2 -1
    h/(2*sina) -1 0 0 0 0 b/2 1
    0 h/sina 0 a 0 h/sina 0 b
    h^3/(12*sina^2) 0 0 a*h/2 -h^3/(12*sina^2) 0 0 -b*h/2
    h^3*cosa/(12*sina^3) (a+b)*h/(2*sina) -a^3/12 0  -h^3*cosa/(12*sina^3) (a+b)*h/(2*sina) b^3/12 0
    1 0 0 0 0 0 0 0];
   
f=[0 0 0 0 0 0 0 1]';
d1=[ 1
    -(h*(a^2 - b^2))/(2*sina*(a^2 + b^2))
    -(2*b^2*h)/(a*sina*(a^2 + b^2))
                                     0
                                     1
  (h*(a^2 - b^2))/(2*sina*(a^2 + b^2))
       -(2*a^2*h)/(b*sina*(a^2 + b^2))
                                     0]
d= k\f