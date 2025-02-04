clc;clear;
syms s a b h t E1 G l t cosa sina real
syms RSs1  RSs2  RSs3  RSs4 RNs1  RNs2  RNs3  RNs4 RZs WSs WNs WZs1 WZs2 WZs3 WZs4 DUs DSs1 DSs2 DSs3 DSs4 DNs1 DNs2 DNs3 DNs4 real
syms dWZs1ds dWZs2ds dWZs3ds dWZs4ds real
syms a0 a1 a2 a3 b0 b1 b2 b3 c0 c1 c2 c3 d0 d1 d2 d3 real

a=75*1e-3;
b=25*1e-3;
h=25*1e-3;
alpha=atan(2*h/(a-b));
angle= alpha*180/pi;
sina=sind(angle);
cosa=cosd(angle);
t=1e-3;
l=750*1e-3;
E=200;
NU=0.3;
E1=E/(1-NU^2);
G = E/(2*(1+NU));

a0= (cosa*h*(4*a^2*b*sina^2 + 11*a^2*h*sina + 4*a*b^2*sina^2 + 10*a*b*h*sina + 12*a*h^2 + 11*b^2*h*sina + 12*b*h^2))/(8*a*b*sina*(3*h^2 + 2*a*h*sina + 2*b*h*sina + a*b*sina^2));
a1=-(cosa*(a - b)*(18*h^2 + 9*a*h*sina + 9*b*h*sina + 4*a*b*sina^2))/(4*a*b*(3*h^2 + 2*a*h*sina + 2*b*h*sina + a*b*sina^2));
a2=-(3*cosa*sina^2*(a - b)^2)/(2*a*b*(3*h^2 + 2*a*h*sina + 2*b*h*sina + a*b*sina^2));
a3=(cosa*sina^2*(a - b)*(6*h + a*sina + b*sina))/(a*b*h*(3*h^2 + 2*a*h*sina + 2*b*h*sina + a*b*sina^2));
b0=1;
b1=(cosa*sina*(3*h + b*sina)*(a - b))/(2*b*(3*h^2 + 2*a*h*sina + 2*b*h*sina + a*b*sina^2));
b2=0;
b3=-(2*cosa*sina*(3*h + b*sina)*(a - b))/(a^2*b*(3*h^2 + 2*a*h*sina + 2*b*h*sina + a*b*sina^2));
c0=-(cosa*h*(4*a^2*b*sina^2 + 11*a^2*h*sina + 4*a*b^2*sina^2 + 10*a*b*h*sina + 12*a*h^2 + 11*b^2*h*sina + 12*b*h^2))/(8*a*b*sina*(3*h^2 + 2*a*h*sina + 2*b*h*sina + a*b*sina^2));
c1=-(cosa*(a - b)*(18*h^2 + 9*a*h*sina + 9*b*h*sina + 4*a*b*sina^2))/(4*a*b*(3*h^2 + 2*a*h*sina + 2*b*h*sina + a*b*sina^2));
c2=(3*cosa*sina^2*(a - b)^2)/(2*a*b*(3*h^2 + 2*a*h*sina + 2*b*h*sina + a*b*sina^2));
c3=(cosa*sina^2*(a - b)*(6*h + a*sina + b*sina))/(a*b*h*(3*h^2 + 2*a*h*sina + 2*b*h*sina + a*b*sina^2));
d0=-1;
d1=(cosa*sina*(3*h + a*sina)*(a - b))/(2*a*(3*h^2 + 2*a*h*sina + 2*b*h*sina + a*b*sina^2));
d2=0;
d3=-(2*cosa*sina*(3*h + a*sina)*(a - b))/(a*b^2*(3*h^2 + 2*a*h*sina + 2*b*h*sina + a*b*sina^2));
DUs=0;
DSs1=1;
DSs2=-h*cosa/(a*sina);
DSs3=-1;
DSs4=h*cosa/(b*sina);
DNs1=a0+a1*s+a2*s^2+a3*s^3;
DNs2=b0+b1*s+b2*s^2+b3*s^3;
DNs3=c0+c1*s+c2*s^2+c3*s^3;
DNs4=d0+d1*s+d2*s^2+d3*s^3;
RSs1=(a+b)*sina/4;
RSs2=h/2;
RSs3=(a+b)*sina/4;
RSs4=h/2;
RNs1=s-(a+b)*cosa/4;
RNs2=s;
RNs3=s-(a+b)*cosa/4;
RNs4=s;
RZs=0;
WSs=0;
WNs=0;
WZs1=s-h*(a^2-b^2)/(2*sina*(a^2+b^2));
WZs2=s*(-2*h*b^2)/(a*sina*(a^2+b^2));
WZs3=s+h*(a^2-b^2)/(2*sina*(a^2+b^2));
WZs4=s*(-2*h*a^2)/(b*sina*(a^2+b^2));
dWZs1ds=1;
dWZs2ds=(-2*h*b^2)/(a*sina*(a^2+b^2));
dWZs3ds=1;
dWZs4ds=(-2*h*a^2)/(b*sina*(a^2+b^2));
d2DNs1ds2=a1+2*a2*s+3*a3*s^2;
d2DNs2ds2=b1+2*b2*s+3*b3*s^2;
d2DNs3ds2=c1+2*c2*s+3*c3*s^2;
d2DNs4ds2=d1+2*d2*s+3*d3*s^2;

A=int(WZs1^2,s,0,h/sina)+int(WZs2^2,s,0,a)+int(WZs3^2,0,h/sina)+int(WZs4^2,s,0,b);
B1star=int(RSs1^2,s,0,h/sina)+int(RSs2^2,s,0,a)+int(RSs3^2,s,0,h/sina)+int(RSs4^2,s,0,b);
B1=int(dWZs1ds^2,s,0,h/sina)+int(dWZs2ds^2,s,0,a)+int(dWZs3ds^2,s,0,h/sina)+int(dWZs4ds^2,s,0,b);
B2=int(RSs1*dWZs1ds,s,0,h/sina)+int(RSs2*dWZs2ds,s,0,a)+int(RSs3*dWZs3ds,s,0,h/sina)+int(RSs4*dWZs4ds,s,0,b);
B3=int(DSs1*dWZs1ds,s,0,h/sina)+int(DSs2*dWZs2ds,s,0,a)+int(DSs3*dWZs3ds,s,0,h/sina)+int(DSs4*dWZs4ds,s,0,b);
B4=int(DSs1*RSs1,s,0,h/sina)+int(DSs2*RSs2,s,0,a)+int(DSs3*RSs3,s,0,h/sina)+int(DSs4*RSs4,s,0,b);
B5=int(DSs1^2,s,0,h/sina)+int(DSs2^2,s,0,a)+int(DSs3^2,s,0,h/sina)+int(DSs4^2,s,0,b);
C=int(d2DNs1ds2^2,s,0,h/sina)+int(d2DNs2ds2^2,s,0,a)+int(d2DNs3ds2^2,s,0,h/sina)+int(d2DNs4ds2^2,s,0,b);

K=[   G*B1star/l      -G*B2/2                 G*B4/l                     -G*B1star/l       -G*B2/2                  -G*B4/l;
          -G*B2/2       (G*B1*l/3+E1*a/l)    -G*B3/2                     G*B2/2      (G*B1*l/6-E1*a/l)      -G*B3/2;
           G*B4/l      -G*B3/2                (G*B5/l+E1*C*l/3)      -G*B4/l      -G*B3/2                    (-G*B5/l+E1*C*l/6);
           -G*B1star/l      G*B2/2                 -G*B4/l                    -G*B1star/l      G*B2/2                     G*B4/l;
           -G*B2/2      (G*B1*l/6-E1*a/l)   -G*B3/2                     G*B2/2       (G*B1*l/3+E1*a/l)       G*B3/2;
           -G*B4/l     -G*B3/2                 (-G*B5/l+E1*C*l/6)       G*B4/l       G*B3/2                     (G*B5/l+E1*C*l/3) ]     
