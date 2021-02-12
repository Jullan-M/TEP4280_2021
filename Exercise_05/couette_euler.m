% Solving the 1D diffusion equation: Start-up of Couette flow 
% using Forward-Time-Central-Space (FTCS) scheme,
% and comparing with analytical solutions.
% The equation is non-dimensional.
clear all
close all
clc

set(groot, 'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

% Number of spatial grid points and cell size
jmax=21;
dy=1/(jmax-1)

% Initialize velocity arrays and set right B.C.
u=zeros(jmax,1);
U3=zeros(1,jmax);

u(jmax)=1;

% y-array for computation of exact solution
y=linspace(0,1,jmax);

% Time step and end of simulation
dt=0.0013
tstop=1; % tstop as found in task a)
nmax=ceil(tstop/dt);

% Scale dt stop at tstop
dt=tstop/nmax;
% Stable for r <= 0.5
r=dt/dy^2
for n=1:nmax
    % Do the explicit euler iteration as an array operation
    u(2:jmax-1) = u(2:jmax-1)*(1-2*r)+(u(3:jmax)+u(1:jmax-2))*r;
    
    t=dt*n;
    inv_t=0.5/sqrt(t);
    i = (0:5)'; % Sum the analytic solution (3) up to n=5
    U3= sum(erfc((2*i+1-y)*inv_t)-erfc((2*i+1+y)*inv_t)); % Analytic solution
    plot(u,y,U3,y,'m+')
    title(['$t=', num2str(t),'$'])
    drawnow;
    pause;
end

hold on
xlabel('$u(y)$')
ylabel('$y$')
% plot the steady state solution
plot(y,y,'r-.')
grid()
legend('FTCS','Analytical (3)','Steady state','Location','best')
hold off