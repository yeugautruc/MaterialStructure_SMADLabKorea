clear; clc;
syms t b h s rn rs13 rs24 uz13 uz24 xn13 xn24 b1 b2 xs13 xs24 real
syms a1 a2 a3 a4 a5 a6 a7 a8 a9 real
syms INSIDE E1 G wa ro di z L N1 N2 N1dz N2dz wadz rodz didz real
% given data
La=0.5;
b=25e-3;
h=50e-3;
t=1e-3;
E=200e9;
V=0.3;
E1=E/(1-V^2);
G=76.9e9;
Ne=100; 
L=La/Ne;
%Element stiffness matrix
rn=-s;
rs13=b/2;
rs24=h/2;
uz13=b*s/2;
uz24=-h*s/2;
uz13ds=b/2;
uz24ds=-h/2;
b1=-4/(h*(b+h));
b2=(2*b+h)/(b+h);
xn13=-b1*s^3+b2*s;
xn24=b1*s^3-b2*s;
xn13ds=-3*b1*s^2+b2;
xn24ds=3*b1*s^2-b2;
xn13ds2=-6*b1*s;
xn24ds2=6*b1*s;
xs13=b*h/(b+h);
xs24=-b*h/(b+h);
a1=t*(2*int(uz13,s,0,h)+2*int(uz24,s,0,b));
a2=(t^3/3)*(2*int(xn13ds2^2,s,0,h)+2*int(xn24ds2^2,s,0,b));
a3=t*(2*int(uz13ds^2,s,0,h)+2*int(uz24ds^2,s,0,b));
a4=t*(2*int(rs13^2,s,0,h)+2*int(rs24^2,s,0,b));
a5=t*(2*int(xs13^2,s,0,h)+2*int(xs24^2,s,0,b));
a6=(t^3/3)*(2*int(xn13ds^2,s,0,h)+2*int(xn24ds^2,s,0,b));
a7=t*(2*int(uz13ds*rs13,s,0,h)+2*int(uz24ds*rs24,s,0,b));
a8=t*(2*int(uz13ds*xs13,s,0,h)+2*int(uz24ds*xs24,s,0,b));
a9=t*(2*int(rs13*xs13,s,0,h)+2*int(rs24*xs24,s,0,b));
N1=1-z/L;
N2=z/L;
N1dz=-1/L;
N2dz=1/L;
ro = [N1 0 0 N2 0 0]';
wa =[0 N1 0 0 N2 0]';
di = [0 0 N1 0 0 N2]';
rodz = [N1dz 0 0 N2dz 0 0]';
wadz = [0 N1dz 0 0 N2dz 0]';
didz = [0 0 N1dz 0 0 N2dz]';
INSIDE= 2*E1*wadz*wadz'*a1 + 2*E1*a2*didz*didz'+G*a3*wa*wa'+G*a4*rodz*rodz'+G*a5*didz*didz'+...
G*a6*di*di'+G*a7*wa*rodz'+G*a7*rodz*wa'+G*a8*wa*di'+G*a8*di*wa'+G*a9*rodz*didz'+G*a9*didz*rodz';
ke= int(INSIDE,0,L)/4

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
Kg;
%calculate dof
Kbig =  mean(diag(Kg))*1e7;
F=zeros(3*(Ne+1),1);
F(3*Ne+2,1)=F(3*Ne+2,1)+Kbig*t;
Kg(1,1)=Kg(1,1)+Kbig;
Kg(2,2)=Kg(2,2)+Kbig;
Kg(3,3)=Kg(3,3)+Kbig;
Kg(3*Ne+2,3*Ne+2)=Kg(3*Ne+2,3*Ne+2)+Kbig;

Dof= Kg\F;

U=zeros(Ne+1,1);
for i=1:(Ne+1)
    U(i)=Dof(3*i-1);
end

length=0:L:La;
plot(length,U/t);
