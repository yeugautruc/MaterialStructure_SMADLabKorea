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
a=0;
c=768*t^3/(b+h);
b1=b*h*(b+h)*t/2;
b2=2*t*(b*h)^2/(h+b);
b3=t^3/3*(276*b*h/(b+h))+134*(b+h);
b4=b*h*(b-h)*t/2;
b5=0;
b6=2*b^2*h^2*t/(h+b);
Ne=20 ; 
Le=L/Ne;
%Number of elements
% element stiffness matrix
Kel=Le*[   G*b1/Le      -G*b4/2                 G*b5/Le                     -G*b1/Le       -G*b4/2                  -G*b5/Le;
          -G*b4/2       (G*b1*Le/3+E*a/Le)    -G*b6/2                     G*b4/2      (G*b1*Le/6-E1*a/Le)      -G*b6/2;
           G*b5/Le      -G*b6/2                 (G*(b2+b3)/Le+E1*c*Le/3)      -G*b5/Le      -G*b6/2                    (-G*(b2+b3)/Le+E1*c*Le/6);
           -G*b1/Le      G*b4/2                 -G*b5/Le                    -G*b1/Le      G*b4/2                     G*b5/Le;
           -G*b4/2      (G*b1*Le/6-E1*a/Le)   -G*b6/2                     G*b4/2      (G*b1*Le/3+E1*a/Le)       G*b6/2;
           -G*b5/Le     -G*b6/2                 -G*(b2+b3)/Le+E1*c*Le/6       G*b5/Le       G*b6/2                     (G*(b2+b3)/Le+E1*c*Le/3) ];
Kg=zeros(3*(Ne+1),3*(Ne+1));
 s=0;
for n=1:Ne
    for i=1:6
        for j=1:6
            Kg(i+s,j+s)= Kg(i+s,j+s)+Kel(i,j);
        end 
    end
    s=s+3;
end
Kg
%calculate dof
Kbig =  mean(diag(Kg))*1e8;
F=zeros(3*(Ne+1),1);

F(3*Ne+1,1)=F(3*Ne+1,1)+100;
Kg(1,1)=Kg(1,1)+Kbig;
Kg(2,2)=Kg(2,2)+Kbig;
Kg(3,3)=Kg(3,3)+Kbig;
%Kg(3*Ne+1,1)=Kg(3*Ne+1,1)+Kbig;
Kg(3*Ne+2,3*Ne+2)=Kg(3*Ne+2,3*Ne+2)+Kbig;
Kg(3*Ne+3,3*Ne+2)=Kg(3*Ne+3,3*Ne+2)+Kbig;
Dof= Kg\F;
U=zeros(Ne+1);
for i=1:(Ne+1)
    U(i)=Dof(3*i-1);
end
Ufind=zeros(Ne+1,1);
for i=1:Ne+1
    Ufind(i)=U(i)*(b*h/4)/t;
end

length=0:Le:0.5;
plot(length,Ufind);