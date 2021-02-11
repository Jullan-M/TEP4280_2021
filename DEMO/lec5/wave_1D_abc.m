% CTCS scheme for 1D wave equation u_tt = c^2 u_xx, where c is speed of sound, 
% with absorbing boundary conditions (ABC), also called nonreflecting BCs.
% The boundaries are artificial, because the acoustic waves are supposed
% to propagate to infinity.


clear all
close all
clc

xmax=10;

jmax=101; 
x=linspace(-xmax,xmax,jmax);
dx=x(2)-x(1);

u=exp(-x.^2);
nu=u;
uo=u;

C=1;
tmax=15;

c=1;
dt=dx*C/c;
nmax=ceil(tmax/dt);
dt=tmax/nmax;
C=c*dt/dx;

plot(x,u)
axis([-xmax xmax -1 1])
pause

j=[2:jmax-1];

tic
for n=1:nmax
    % CTCS
    un(j)=2*u(j)-uo(j)+C^2*(u(j+1)-2*u(j)+u(j-1));
    
    % Nonreflecting boundary conditions
    un(1)=(1-C)*u(1)+C*u(2);
    un(jmax)=(1-C)*u(jmax)+C*u(jmax-1);
    
    uo=u;   % set new u^{n-1} equal old u^{n}
    u=un;   % set new u^{n} equal old u^{n+1}

    plot(x,u)
    axis([-xmax xmax -1 1])
    title(['Time=',num2str(n*dt,2)])
    drawnow;
end
toc

xlabel('x')
ylabel('u')