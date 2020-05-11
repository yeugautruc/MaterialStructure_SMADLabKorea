clc;
clear;
close;
K=44.5;
L=3;
M=3;
h=-1000;
g=293.15
fc=1000;
n=2
N=n^2;
Le=L/n;
Me=M/n;
elnode=[];
for i=1:N
    p=floor(Le*(i-1)/L);
    elnode=[elnode;[i+p i+p+1 i+p+1+n+1 i+p+1+n]];
end
ek=[K*Le/(3*Me)+K*Me/(3*Le) K*Le/(6*Me)-K*Me/(3*Le) -K*Le/(6*Me)-K*Me/(6*Le) K*Le/(6*Me)-K*Me/(3*Le)
0 K*Le/(3*Me)+K*Me/(3*Le) K*Le/(6*Me)-K*Me/(3*Le) -K*Le/(6*Me)-K*Me/(6*Le)
0 0 K*Le/(3*Me)+K*Me/(3*Le) K*Le/(6*Me)-K*Me/(3*Le)
0 0 0 K*Le/(3*Me)+K*Me/(3*Le)];
ek=ek+ek';
for i=1:length(ek)
    ek(i,i)=ek(i,i)/2;
end


fb=[Le*Me/4
    Le*Me/4
    Me^2/4
    Me^2/4]*fc;
fe=[0
    0
    Me/2
    Me/2]*h;
k=zeros((n+1)^2);
f=zeros((n+1)^2,1);
for i=1:N
    p=floor(Le*(i-1)/L);
    dof=elnode(i,:);
    k(dof,dof)=k(dof,dof)+ek;
    f(dof)=f(dof)+fb;
    if p==n-1
        f(dof)=f(dof)-fe;
    end
end
fixdof=[1:n+1:(n+1)^2];
kbig=max(diag(k))*1e7;
for i=fixdof
    k(i,i)=k(i,i)+kbig;
    f(i)=f(i)+g*kbig;
end
d=k\f

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