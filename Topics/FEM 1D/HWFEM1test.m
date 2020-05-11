clear all
clc

Ne=10;
le=1/Ne;
eqn = 'D2u + x^2 = 0';
inits = 'Du(0) = 2, u(1) = 1';
x = 0:1/Ne:1-1/Ne;
y = dsolve(eqn,inits,'x')
z = eval(vectorize(y));
Ke = 1/le*[1 -1;
        -1 1];
Kg=zeros(Ne,Ne);
 s=0;
for n=1:Ne-1
    for i=1:2
        for j=1:2
            Kg(i+s,j+s)= Kg(i+s,j+s)+Ke(i,j);
        end 
    end
    s=s+1;
end
%Fe=le^3*[ 1/12;1/4];
Fg=zeros(Ne,1);

 for i=1: Ne-1
     %Fe1=1/(4*le);
     %Fe2=(1/3-1/(4*le));
     Fe1=(((i+1)^4-i^4-4*i^3)/(12*Ne^3));
     Fe2=(((i+1)^4-i^4)/(4*Ne^3));
     Fg(i,1)=Fg(i,1)+Fe1;
     Fg(i+1,1)= Fg(i+1)+Fe2;
 end

%Kbig =  mean(diag(Kg))*1e5;
%Kg(1,1)=Kg(1,1)+Kbig;
%Kg(Ne+1,Ne+1)=Kg(Ne+1,Ne+1)+Kbig;
Fg(1,1)=Fg(1,1)-2;
Fg(Ne-1,1)=Fg(Ne-1,1)+1/le;
Kg
Fg
d= Kg\Fg;

plot(x,d,x,z)