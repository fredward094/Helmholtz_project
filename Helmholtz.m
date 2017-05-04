clc
clear
% The constants
n=9;
ax=-pi; ay=-pi; bx=pi; by=pi;
gamma=-pi;
% Generating x,y
x=linspace(ax,bx,n); y=linspace(ay,by,n);
% Boundary conditions
phiab=cos(pi*(y-ay)).*cosh(by-y);
psiab=(y-ay).^2.*sin(pi*(y-ay)/(by-ay));
u(:,1)=phiab; u(:,n)=psiab;
h=lenght(x);

F = cos(pi/2*(2*(x-ax)./(bx-ax)+1)).*sin(pi*(y-ay)./(by-ay));

for j=2:n-1
    for i=2:n-1
        F(i,j) = cos(pi/2*(2*(x(i)-ax)./(bx-ax)+1)).*sin(pi*(y(j)-ay)./(by-ay));
        u(i,j)= 1/(4-gamma)*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1))*(1/h);
    end
end
u(1,:)=u(2,2); u(n,:)=u(n-1,n-1);
U=u
mesh(U)
[Solution,Error_estimate,Number_of_iterations,flag]=SOR_trial1(U,zeros(n,1),F1,1,1e4,1e-4)
