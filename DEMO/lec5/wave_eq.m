% CTCS for 1D wave equation
clear all
close all
clc

xmax=10;
jmax=101;
tmax=80;

x=linspace(-xmax,xmax,jmax);

dx= x(2)-x(1); % Grid spacing
u=exp(-x.^2); % % u: u^n Gaussian distribution
u(1)=0; u(jmax)=0; % Dirichlet BCs
C=1; % Courant number
c=1; % Speed of medium (sound)
dt=C*dx/c; % Time step

nmax=ceil(tmax/dt);
uo=u; % uo: u^n-1
un=u; % um: u^n+1

j=[2:jmax-1];
tic
for n=1:nmax
    un(j)=2*u(j)-uo(j)+C^2*(u(j+1)-2*u(j)+u(j-1));
    % Neumann BCs
    un(1)=un(2);
    un(jmax)=un(jmax-1);
    
    uo=u;
    u=un;
    t=dt*n;
    
    plot(x,u,'b')
    axis([-xmax, xmax, -1, 1])
    title(['t = ',num2str(t,2)])
    drawnow;
end
toc

xlabel('x')
ylabel('u')