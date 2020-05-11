clear;
clc;
% given data
L=0.5;
b=25e-3;
h=50e-3;
t=1e-3;
E=200e9;
V=0.3;
E1=E/(1-V^2);
G=76.9e9;
Ne=100; 
Le=L/Ne;
%Element stiffness matrix
a1=(t*b^2*h)/2 + (t*b*h^2)/2;
a2=(b*t^3)/6 + (h*t^3)/6;
a3=(b^2*h*t)/2 - (b*h^2*t)/2;
a4=0;
a5=0;
a6=(t*b^2*h)/2 + (t*b*h^2)/2;
a7=(t*b^3*h^2)/24 + (t*b^2*h^3)/24;
a8=(2*b^2*h^2*t)/(b + h);
a9=(b^2*h*t)/2 - (b*h^2*t)/2;
a10=0;
a11=(2*b^2*h^3*t)/(b + h)^2 + (2*b^3*h^2*t)/(b + h)^2;
a12=(2*b*t^3*(b^2 + 5*h^2))/(15*(b + h)^2) + (2*h*t^3*(5*b^2 + h^2))/(15*(b + h)^2);
a13=(8*b*t^3)/(b + h)^2 + (8*h*t^3)/(b + h)^2;
ke=zeros(6);
ke(1,1)=G/Le*(a1+a2);
ke(1,2)=-G*a3/2;
ke(1,3)=G/Le*(a4-a5);
ke(1,4)=-G/Le*(a1+a2);
ke(1,5)=-G*a3/2;
ke(1,6)=G/Le*(a5-a4);
ke(2,2)=Le*G*a6/3+E1*a7/Le;
ke(2,3)=-G*a8/2;
ke(2,4)=G*a9-G*a9/2;
ke(2,5)=G*Le*a6/6-E1*a7/Le;
ke(2,6)=G*a8/2;
ke(3,3)=G/Le*a11+G/Le*a12+E1*Le*a13/3;
ke(3,4)=G/Le*(a5-a4);
ke(3,5)=-G*a8/2;
ke(3,6)=-G/Le*a11-G/Le*a12+E1*Le*a13/6;
ke(4,4)=G/Le*(a1+a2);
ke(4,5)=G*a3/2;
ke(4,6)=G/Le*(a4-a5);
ke(5,5)=E1/Le*a7+Le*G*a6/3;
ke(5,6)=G*a8/2;
ke(6,6)=G/Le*a11+G/Le*a12+E1*Le*a13/3;
for i=2:6
    for j=1:i-1
        ke(i,j)=ke(j,i);
    end
end
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
    U(i)=Dof(3*i);
end

length=0:Le:L;
plot(length,U/t);
