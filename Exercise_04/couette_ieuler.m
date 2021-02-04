% Solving the 1D diffusion equation: Start-up of Couette flow 
% using Forward-Time-Central-Space (FTCS) scheme (Implicit Euler),
% and comparing with analytical solutions.

% The equation is non-dimensional.
clear all
close all
clc

% Number of grid points (excluding boundary points) and cell size
jmax=21;
dy=1/(jmax-1);

% Initialize velocity arrays and set right B.C.
u=zeros(jmax,1);
U3=zeros(1,jmax);
ua=0; ub=1;

% y-array for computation of exact solution
y=linspace(0,1,jmax);

% Time step and end of simulation
dt=0.00125;
tstop=1;
nmax=ceil(tstop/dt);

% Scale dt stop at tstop
dt=tstop/nmax;
% Stable for r <= 0.5
r=dt/dy^2;
A = diag((1+2*r)*ones(1,jmax)) + diag(-r*ones(1,jmax-1),1) + diag(-r*ones(1,jmax-1),-1);
b = u;

for n=1:nmax
    b(1) = u(1)+r*ua;
    b(end) = u(end)+r*ub;
    u = A \ b;
    b = u;
    
    t=dt*n;
    t=0.5/sqrt(t);
    
    i = (0:5)';
    U3= sum(erfc((2*i+1-y)*t)-erfc((2*i+1+y)*t));
    
    plot(u,y,U3,y,'m+')
    drawnow;
    % pause                   
end

hold on
xlabel('u(y)')
ylabel('y')
% plot the steady state solution
plot(y,y,'r-.')                  
legend('FTCS','Analytical (3)','Steady state','Location','best')
hold off