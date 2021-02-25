% Explicit upwind method for 1D linear advection equation
clear all
close all
clc

set(groot, 'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

% Number of spatial grid points and cell size
jmax=101;
dx=1/(jmax-1);
x=linspace(0,1,jmax);

% Initialize velocity arrays and set right ICs and BC.
k=2*pi;
u0=zeros(1,jmax);
u0(jmax)=1000; % Right boundary at T=1000C at all times
u=u0; unew=u0; u_ex=u0;

% x-array for computation of exact solution
x=linspace(0,1,jmax);

c=-0.8; % Advection velocity
C=-0.5; % Courant number
% Time step and end of simulation
dt=C*dx/c
tmax=1; % Max time step
nmax=ceil(tmax/dt);

% Scale dt stop at tstop
dt=tmax/nmax;
j =[1:jmax-1];
figure(1)
hvsd = @(x) x >= 0; % Modified "Heaviside" function

for n=1:nmax
    unew(j) = u(j)*(1+C)-C*u(j+1);
    u=unew;
    t=n*dt;
    % Analytical solution: Heaviside function propagating to the left
    u_ex = 1000*hvsd(x-1-c*t); 
    
    plot(x, u, 'x', x, u_ex, 'r')
    title(['$t=', num2str(t),'$ s'])
    pause(0.025)
    drawnow;
end

hold on
grid()
xlabel('$x$ [m]')
ylabel('$T$ [$^\circ$ C]')
legend('Explicit Upwind Method','Analytical','Location','best')
hold off