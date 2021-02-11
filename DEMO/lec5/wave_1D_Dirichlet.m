% CTCS scheme for 1D wave equation u_tt = c^2 u_xx, where c is speed of sound, 
% Dirichlet boundary conditions.
% The boundaries of the domain, e.g. of an organ pipe, are supposed to be
% open to the ambient such that the boundary conditions are
% homogeneous Dirichlet boundary conditions, i.e., the acoustic pressure
% u is zero at the boundaries.

clear all
close all
clc

xmax=10;
jmax=101;
x=linspace(-xmax,xmax,jmax);
dx=x(2)-x(1);

u=exp(-x.^2); % initial condition: Gauss pulse

u(1)=0; u(jmax)=0; % Dirichlet BCs

tmax=80;

C=1; % Courant number c dt/dx
c=1;
dt=C*dx/c;
nmax=ceil(tmax/dt);
dt=tmax/nmax;
C=c*dt/dx;

plot(x,u)
axis([-xmax xmax -1 1])
pause

% u: u^n
un=u; % un: u^n+1
uo=u; % uo: u^n-1

j=[2:jmax-1];

tic
for n=1:nmax
    % CTCS
    % using for loop 
%     for j=2:jmax-1
%         un(j)=2*u(j)-uo(j)+C^2*(u(j+1)-2*u(j)+u(j-1));
%     end

    % using vector operation which is faster than for loop in MATLAB  
    un(j)=2*u(j)-uo(j)+C^2*(u(j+1)-2*u(j)+u(j-1));

    uo=u; % u_n-1 is set to u_n
    u=un; % u_n is set to u_n+1
    
    t=n*dt;
    plot(x,u)
    axis([-xmax xmax -1 1])
    title(['Time=',num2str(t,2)])
    drawnow
end
toc

xlabel('x')
ylabel('u')

    
    
    