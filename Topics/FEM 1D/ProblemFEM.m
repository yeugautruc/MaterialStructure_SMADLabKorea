clear all
clc

Ne=4;
le=1/Ne;
%exactly solution
x = 0:1/Ne:1;
xe=0:1/100:1
for i=1:size(xe,2)
y(i) = - xe(i)^4/12 + 2*xe(i) - 11/12;
end
%FEM Solution
%Stiffness matrix
Ke = 1/le*[1 -1;
        -1 1];
%Global matrix
Kg=zeros(Ne+1,Ne+1);
 s=0;
for n=1:Ne
    for i=1:2
        for j=1:2
            Kg(i+s,j+s)= Kg(i+s,j+s)+Ke(i,j);
        end 
    end
   s=s+1;
end

Fg=zeros(Ne+1,1);
 for i=1: Ne
   % Fe1=((i^4-(i-1)^4-4*(i-1)^3)/(12*Ne^3));
   % Fe2=((i^4-(i-1)^4)/(4*Ne^3));
  Fe1=1/le*(le*i/3*((le*i)^3-(le*(i-1))^3)-1/4*((le*i)^4-(le*(i-1))^4));
  Fe2=1/le*(1/4*((le*i)^4-(le*(i-1))^4)-le*(i-1)/3*((i*le)^3-(le*(i-1))^3));
     Fg(i,1)=Fg(i,1)+Fe1;
     Fg(i+1,1)= Fg(i+1)+Fe2;
 end
Fg(Ne)=Fg(Ne)-Kg(Ne,Ne+1)
Fg(1,1)=Fg(1,1)-2;
Kg(Ne+1,:)=[];
Kg(:,Ne+1)=[];
Fg(Ne+1)=[];
Fg
d= Kg\Fg;
d=[d;1];
plot(x,d,'r*--',xe,y,'b-')