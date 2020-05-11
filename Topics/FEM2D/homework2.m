clear;
clc;
%Given data
K=44.5;
f=1000;
h=1000;
m=3/2;
l=3/2;
g=293.15;
%Element stiffness matrix
ek(1,1)=K*(m^2+l^2)/(3*m*l);
ek(1,2)=K*(-2*m^2+l^2)/(6*m*l);
ek(1,3)=K*(-m^2-l^2)/(6*m*l);
ek(1,4)=K*(m^2-2*l^2)/(6*m*l);
ek(2,2)=K*(m^2+l^2)/(3*m*l);
ek(2,3)=K*(m^2-2*l^2)/(6*m*l);
ek(2,4)=K*(-m^2-l^2)/(6*m*l);
ek(3,3)=K*(m^2+l^2)/(3*m*l);
ek(3,4)=K*(-2*m^2+l^2)/(6*m*l);
ek(4,4)=K*(m^2+l^2)/(3*m*l);
for i=2:4
    for j=1:i-1
        ek(i,j)=ek(j,i);
    end
end
ek;
%Etable of mesh
etable=[1 2 5 4;2 3 6 5;4 5 8 7; 5 6 9 8];
%Global stiffness matrix
kglobal=zeros(9);
for i=1:4
  kglobal(etable(i,:),etable(i,:))=kglobal(etable(i,:),etable(i,:))+ek;
end
kglobal
%Element load matrix
fe=f*m*l/4*[1;1;1;1];
%Global load matrix
fglobal=zeros(9,1);
for i=1:4
 %for i=1:9
fglobal(etable(i,:))=fglobal(etable(i,:))+fe;
end

%Apply the BC condition
big=mean(diag(kglobal))*1e6;
fixdof=[1 4 7];
for i=fixdof
    kglobal(i,i)=kglobal(i,i)+big;
    fglobal(i,1)=fglobal(i,1)+g*big;
end
fglobal(7,1)=fglobal(7,1)+l/2*h;
fglobal(8,1)=fglobal(8,1)+l/2*h+l/2*h;
fglobal(9,1)=fglobal(9,1)+l/2*h;
fglobal
%Calculate the solution vector
d=kglobal\fglobal
for i=1:3
    for j=1:3
        u1(i,j)=d(j+3*(i-1),1);
    end 
end
    for i=1:3
        for j=1:3
        u(i,j)=u1(4-i,j);
    end 
    end
u
surf(u)
       
