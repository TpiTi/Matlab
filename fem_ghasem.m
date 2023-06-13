clc;
clear all;
close all;
%% 

E=70e9;
v=0.3;
D=(E/((1+v)*(1-2*v)))*[1-v v v 0;v 1-v v 0;v v 1-v 0;0 0 0 (1-2*v)/2];
r=[1 2 1];
z=[0 1.5 3];
zbar=sum(z)/3;
rbar=sum(r)/3;
A=[1 r(1) z(1);1 r(2) z(2);1 r(3) z(3)];
Aa=A(:,2:3);
cros_section=det(A)/2;
alpha=zeros(1,3);
beta=[(z(2)-z(3)) (z(3)-z(1)) (z(1)-z(2))];
gama=[(r(3)-r(2)) (r(1)-r(3)) (r(2)-r(1))];
for i=1:3
    m=Aa;
    m(i,:)=[];
alpha(i)=((-1)^(i+1))*det(m);

B(1:4,2*i-1:2*i)=(1/det(A))*[beta(i) 0;0 gama(i);(alpha(i)/rbar)+beta(i)+(gama(i)*zbar/rbar) 0;gama(i) beta(i)]; 
end


for i=1:3
   K{i}=(2*pi)*(rbar*zbar)*(B(:,2*i-1:2*i))'*D*(B(:,2*i-1:2*i));
end
g=zeros(2);
Ktotal=[K{1} g g;g K{2} g;g g K{3}];
q=[47.1;0;0;0;47.1;0];
x=inv(Ktotal)*q;
strain=(1/det(A))*B*x;
stress=D*strain;





