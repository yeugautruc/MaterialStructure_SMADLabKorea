clc;clear;
syms s a b h t E1 G l t cosa sina real
syms RSs1  RSs2  RSs3  RSs4 RNs1  RNs2  RNs3  RNs4 RZs WSs WNs WZs1 WZs2 WZs3 WZs4 DUs DSs1 DSs2 DSs3 DSs4 DNs1 DNs2 DNs3 DNs4 real
syms dWZs1ds dWZs2ds dWZs3ds dWZs4ds real
syms a0 a1 a2 a3 b0 b1 b2 b3 c0 c1 c2 c3 d0 d1 d2 d3 real
syms A B1star B1 B2 B3 B4 B5 C z real
% given data
La=0.5;
Ne=50; 
l=La/Ne;
a=75*1e-3;
b=25*1e-3;
h=25*1e-3;
alpha=atan(2*h/(a-b));
angle= alpha*180/pi;
sina=sind(angle);
cosa=cosd(angle);
t=1e-3;
E=200;
NU=0.3;
E1=E/(1-NU^2);
G = E/(2*(1+NU));
%Element stiffness matrix
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
DNs1ds=2*a2+6*a3*s;
DNs2ds=2*b2+6*b3*s;
DNs3ds=2*c2+6*c3*s;
DNs4ds=2*d2+6*d3*s;
d2DNs1ds2=a1+2*a2*s+3*a3*s^2;
d2DNs2ds2=b1+2*b2*s+3*b3*s^2;
d2DNs3ds2=c1+2*c2*s+3*c3*s^2;
d2DNs4ds2=d1+2*d2*s+3*d3*s^2;

A1=t*(int(WZs1,s,0,h/sina)+int(WZs2,s,0,a)+int(WZs3,s,0,h/sina)+int(WZs4,s,0,b));
A2=(t^3/3)*(int(d2DNs1ds2^2,s,0,h/sina)+int(d2DNs2ds2^2,s,0,a)+int(d2DNs3ds2^2,s,0,h/sina)+int(d2DNs4ds2^2,s,0,a));
A3=t*(int(dWZs1ds^2,s,0,h/sina)+int(dWZs2ds^2,s,0,a)+int(dWZs3ds^2,s,0,h/sina)+int(dWZs4ds^2,s,0,b));
A4=t*(int(RSs1^2,s,0,h/sina)+int(RSs2^2,s,0,a)+int(RSs3^2,s,0,h/sina)+int(RSs4^2,s,0,b));
A5=t*(int(DSs1^2,s,0,h/sina)+int(DSs2^2,s,0,a)+int(DSs3^2,s,0,h/sina)+int(DSs4^2,s,0,b));
A6=(t^3/3)*(int(DNs1ds^2,s,0,h/sina)+int(DNs2ds^2,s,0,a)+int(DNs3ds^2,s,0,h/sina)+int(DNs3ds^2,s,0,b));
A7=t*(int(dWZs1ds*RSs1,s,0,h/sina)+int(dWZs2ds*RSs2,s,0,a)+int(dWZs3ds*RSs1,s,0,h/sina)+int(dWZs4ds*RSs2,s,0,b));
A8=t*(int(dWZs1ds*DSs1,s,0,h/sina)+int(dWZs2ds*DSs2,s,0,a)+int(dWZs3ds*DSs3,s,0,h/sina)+int(dWZs4ds*DSs4,s,0,b));
A9=t*(int(RSs1*DSs1,s,0,h/sina)+int(RSs2*DSs2,s,0,a)+int(RSs3*DSs3,s,0,h/sina)+int(RSs4*DSs4,s,0,b));
N1=1-z/l;
N2=z/l;
N1dz=-1/l;
N2dz=1/l;
ro = [N1 0 0 N2 0 0]';
wa =[0 N1 0 0 N2 0]';
di = [0 0 N1 0 0 N2]';
rodz = [N1dz 0 0 N2dz 0 0]';
wadz = [0 N1dz 0 0 N2dz 0]';
didz = [0 0 N1dz 0 0 N2dz]';
INSIDE= 2*E1*wadz*wadz'*A1 + 2*E1*A2*didz*didz'+G*A3*wa*wa'+G*A4*rodz*rodz'+G*A5*didz*didz'+...
G*A6*di*di'+G*A7*wa*rodz'+G*A7*rodz*wa'+G*A8*wa*di'+G*A8*di*wa'+G*A9*rodz*didz'+G*A9*didz*rodz';
ke= int(INSIDE,z,0,l)/4;
   %Global stiffness matrix
Kg=zeros(3*(Ne+1),3*(Ne+1));
 s=0;
for n=1:Ne
    for i=1:6
        for j=1:6
            Kg(i+s,j+s)= Kg(i+s,j+s)+ke(i,j);
        end 
    end
    s=s+3;
end
%calculate dof
Kbig =  mean(diag(Kg))*1e8;
F=zeros(3*(Ne+1),1);

F(3*Ne+1,1)=F(3*Ne+1,1)+100;
Kg(1,1)=Kg(1,1)+Kbig;
Kg(2,2)=Kg(2,2)+Kbig;
Kg(3,3)=Kg(3,3)+Kbig;
Kg(3*Ne+2,3*Ne+2)=Kg(3*Ne+2,3*Ne+2)+Kbig;
Kg(3*Ne+3,3*Ne+2)=Kg(3*Ne+3,3*Ne+2)+Kbig;
Dof= Kg\F;

U=zeros(Ne+1,1);
for i=1:(Ne+1)
    U(i)=Dof(3*i-1);
end
Umax=U(1);
for i=1:(Ne+1)
    if U(i)<Umax
        Umax=U(i);
    end
end
length=0:l:La;
plot(length,U/Umax);
