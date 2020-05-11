clear;
clc;
E=200*1e9;
Nu=0.3;
b=50/1000;
h=50/1000;
L=1;
A=b*h;
G=(E/(2*(1+Nu)));
i=b*(h^3)/12;
%Vmax=zeros(10,1);
%en=zeros(100,1);
%en(1)=1;

m= 3
    
%m=8;
eL=L/m;
%cal stiffness
 K0=[1/eL 1/2 -1/eL 1/2;1/2 (eL/3)+((E*i)/(G*A*eL)) -1/2 (eL/6)-((E*i)/(G*A*eL));
-1/eL -1/2 1/eL -1/2;1/2 (eL/6)-((E*i)/(G*A*eL)) -1/2 (eL/3)+((E*i)/(G*A*eL))];
Kt=G*A*K0;
P=1e4;
KG=zeros(2*m+2,2*m+2);
s=0;
for n=1:m
    for i=1:4
        for j=1:4
            KG(i+s,j+s)= KG(i+s,j+s)+Kt(i,j);
        end 
    end
    s=s+2;
end
%Global stiff matrix
Force = zeros(2*m+2,1);
Force(2*m+1,1)= P;
Pen = mean(diag(KG))*1e78; %Gauss
KG(1,1) = KG(1,1) +Pen;
KG(2,2) = KG(2,2) +Pen;
Dist= KG\Force;
V=zeros(m+1,1);
for i=1:m+1
    V(i,1)=Dist(2*i-1,1);
end
%Vmax(m,1)=V(m+1,1)
%cal displacement
%length=zeros(m+1,1);
%for i=1:m+1
   % length(i,1)=(L/m)*(i-1);
%end
%V=zeros(m+1,1);
%for i=1:79
 %  en(i+1,1)=en(i,1)+1;
%end
%end
%Vmax(t)

%end
%plot(length,V);
%en;
%plot(en,Vmax);





    