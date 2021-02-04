% Solving the 1D diffusion equation with source term: 
% Start-up of Poiseulle flow 
% using Forward-Time-Central-Space (FTCS) scheme,
% and comparing with analytical steady state solution.
% The equation is non-dimensional.

clear all
close all
clc

% Number of grid points and cell size
jmax=21;
dy=2/(jmax-1);

% Initialize velocity array and set B.C.s
u=zeros(jmax,1);
U=0; % For U>0, the startup of a combined Poiseuille-Couette flow can be calculated.
u(jmax)=U;

% y-array for computation of exact solution
y=linspace(-1,1,jmax);

% End time of simulation
tstop=3;

% Pressure gradient coefficient k=-(h^3/(rho*nu^2))*dp/dx,
% where h is the channel half height
k=2;

% Time step size and number of time steps
% FTCS stable for r <= 0.5
r=0.5;
dt=r*dy^2;
nmax=ceil(tstop/dt);
dt=tstop/nmax;
r=dt/dy^2;

j=[2:jmax-1];

tic
for n=1:nmax
    % Vector operations in MATLAB execute much faster than for loops! 
    u(j)=k*dt + u(j)*(1-2*r)+(u(j+1)+u(j-1))*r; % j is a vector here!

    t=dt*n;
    plot(u,y,'bx')
    xlim([0 1+U])
    title(['t = ',num2str(t,2)])
    drawnow;
    
    % pause         % Use this if you want stepwise simulation                          
end
toc

hold on
xlabel('u(y)')
ylabel('y')
% plot the steady state solution
% u_ex=1-y.^2; % for Poiseuille flow
u_ex=0.5*(U+k)+0.5*U*y-0.5*k*y.^2; % for Poiseuille-Couette flow
plot(u_ex,y,'r-.')                  
legend('FTCS','exact steady state solution','Location','best')
hold off

disp(['||u - u_exact steady state||_2=',num2str(norm(u-u_ex')*sqrt(dy))])