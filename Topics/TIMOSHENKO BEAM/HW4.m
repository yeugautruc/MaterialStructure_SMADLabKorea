clear;
clc;
%% Given
E=200e9;
Nu=0.3;
b=20e-3;
h=50e-3;
A=b*h;
G=(E/(2*(1+Nu)));
im=b*(h^3)/12;
L=1;
P=-1000;
Ne=20;
m=Ne;
eL=L/m;
Force = zeros(2*m+2,1);
Force(2*m+1,1)= P;
%% Stiffness Euler
Ke=[12 6*eL -12 6*eL;
    6*eL 4*eL*eL -6*eL 2*eL*eL; 
    -12 -6*eL 12 -6*eL;
    6*eL 2*eL*eL -6*eL 4*eL*eL];
Kee=((E*im)/(eL^3))*Ke;
Kge=zeros(2*m+2,2*m+2);
s=0;
for n=1:m
    for i=1:4
        for j=1:4
            Kge(i+s,j+s)= Kge(i+s,j+s)+Kee(i,j);
        end 
    end
    s=s+2;
end

%% Stiffness timoshenko
K0=[1/eL 1/2 -1/eL 1/2;1/2 (eL/3)+((E*im)/(G*A*eL)) -1/2 (eL/6)-((E*im)/(G*A*eL));
-1/eL -1/2 1/eL -1/2;1/2 (eL/6)-((E*im)/(G*A*eL)) -1/2 (eL/3)+((E*im)/(G*A*eL))];
Kt=G*A*K0;
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
%% Cal for euler
big1 = mean(diag(Kge))*1e7; %Gauss
Kge(1,1) = Kge(1,1) +big1;
Kge(2,2) = Kge(2,2) +big1;

Diste= Kge\Force;
Ve=zeros(m+1,1);
for i=1:m+1
    Ve(i,1)=Diste(2*i-1,1);
end
Vmaxe=max(Ve)

%% Cal for timoshenko
big2 = mean(diag(KG))*1e7; %Gauss
KG(1,1) = KG(1,1) +big2;
KG(2,2) = KG(2,2) +big2;
Distt= KG\Force;
Vt=zeros(m+1,1);
for i=1:m+1
    Vt(i,1)=Distt(2*i-1,1);
end
Vmaxt=max(Vt)
%% Post
length=0:L/Ne:L;
plot(length,Ve,'*-r',length,Vt,'o-b');legend('Euler-Bernouli','Timoshenko');
xlabel('Length');ylabel('Displacement');






    