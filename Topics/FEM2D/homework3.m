clc;
clear;
E=200e9;
NU=0.3;
thk=5e-3;
f=10e6;
R=20e-3;
%Node coodinates
node=[0    0
          0.6 0
          2    0
          0    1
          1    1
          2    1
          0    2
         1.4  2
          2    2];
%Element connectivity
etable=[1 2 5 4
            2 3 6 5
            4 5 8 7
            5 6 9 8];
%Element_dof connectivity
etable_dof=[2*etable(:,1)-1, 2*etable(:,1), 2*etable(:,2)-1, 2*etable(:,2),2*etable(:,3)-1, 2*etable(:,3), 2*etable(:,4)-1, 2*etable(:,4)];

D=E/(1-NU^2)*[1 NU 0;NU 1 0;0 0 (1-NU)/2];
K=sparse(18,18);
F=sparse(18,1);
gp=[-1/sqrt(3), 1/sqrt(3)];
wt=[1,1];
EN=size(etable,1);

for i=1:EN
    %element stiffness matrix calculation
ek=zeros(8);
x=node(etable(i,:),1);
y=node(etable(i,:),2);
for j=1:2
    for k=1:2
    ksi=gp(k);
    eta=gp(j);
    dndksi=[-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4];
    dndeta=[-(1-ksi)/4, -(1+ksi)/4, (1+ksi)/4, (1-ksi)/4];
    dxdksi=dndksi*x;
    dxdeta=dndeta*x;
    dydksi=dndksi*y;
    dydeta=dndeta*y;
    jacob=dxdksi*dydeta-dydksi*dxdeta;
    dn1dx=1/jacob*(dydeta*dndksi(1)-dydksi*dndeta(1));
    dn2dx=1/jacob*(dydeta*dndksi(2)-dydksi*dndeta(2));
    dn3dx=1/jacob*(dydeta*dndksi(3)-dydksi*dndeta(3));
    dn4dx=1/jacob*(dydeta*dndksi(4)-dydksi*dndeta(4));
    dn1dy=1/jacob*(-dxdeta*dndksi(1)+dxdksi*dndeta(1));
    dn2dy=1/jacob*(-dxdeta*dndksi(2)+dxdksi*dndeta(2));
    dn3dy=1/jacob*(-dxdeta*dndksi(3)+dxdksi*dndeta(3));
    dn4dy=1/jacob*(-dxdeta*dndksi(4)+dxdksi*dndeta(4));
    B=[dn1dx 0 dn2dx 0 dn3dx 0 dn4dx 0;
        0 dn1dy 0 dn2dy 0 dn3dy 0 dn4dy;
        dn1dy dn1dx dn2dy dn2dx dn3dy dn3dx dn4dy dn4dx];
       ek=ek+B'*D*B*jacob*thk*wt(j)*wt(k);   
    end 
end
%assemble
edof=etable_dof(i,:);
K(edof,edof)=K(edof,edof)+ek;
end

for i=1:EN
  %force element calculation
fe=zeros(8,1);
x=node(etable(i,:),1);
y=node(etable(i,:),2);
for j=1:2
    for k=1:2
    ksi=gp(k);
    eta=gp(j);
    dndksi=[-(1-eta)/4, (1-eta)/4, (1+eta)/4, -(1+eta)/4];
    dndeta=[-(1-ksi)/4, -(1+ksi)/4, (1+ksi)/4, (1-ksi)/4];
    dxdksi=dndksi*x;
    dxdeta=dndeta*x;
    dydksi=dndksi*y;
    dydeta=dndeta*y;
    jacob=dxdksi*dydeta-dydksi*dxdeta;
       fe=fe+f*jacob*thk*wt(j)*wt(k);   
    end 
end
%assemble
edof=etable_dof(i,:);
F(edof)=F(edof)+fe;
end
spy(K)

