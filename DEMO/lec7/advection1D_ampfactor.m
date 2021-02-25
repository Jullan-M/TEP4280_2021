% Different methods for 1D linear advection equation with periodic boundary conditions
% using their amplification factors
clear all
close all
clc

set(groot, 'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

% Number of spatial grid points and cell size
jmax=21;
dx=1/(jmax-1);
x=linspace(0,1,jmax);

% Initialize velocity arrays and set right B.C.
k=4*pi;
u0=cos(k*x);
u=u0; g=u0; u_ex=u0;

% y-array for computation of exact solution
x=linspace(0,1,jmax);

c=1; % Advection velocity
C=0.25; % Courant number
% Time step and end of simulation
dt=C*dx/c
tmax=0.5; % tstop as found in task a)
nmax=ceil(tmax/dt);

% Scale dt stop at tstop
dt=tmax/nmax;

% Schemes: 
% 1=FTCS, 2=Expl.Upwind, 3=Impl.Upwind, 4=Lax-Friedrichs, 5=Lax-Wendroff
scheme=5;

figure(1)
% Time loop 
for n=1:nmax
    switch scheme
        case 1
            % FTCS
            g=1-i*C*sin(k*dx);
        case 2
            % Explicit upwind method
            g=1-C*(1-exp(-i*k*dx));
        case 3
            % Implicit upwind method
            g=1/(1+C*(1-exp(-i*k*dx)));
        case 4
            % Lax-Friedrichs method
            g=cos(k*dx)-i*C*sin(k*dx);
        case 5
            % Lax-Wendroff method
            g=1-C^2*(1-cos(k*dx))-i*C*sin(k*dx);
    end
    u=real(g^n*exp(i*k*x));
    t=n*dt;
    u_ex = cos(k*(x-c*t));
    
    plot(x, u, 'x', x, u_ex, 'r')
    title(['$t=', num2str(t),'$', 'Scheme = ', num2str(scheme)])
    
    drawnow;
end

hold on
grid()
xlabel('$x$')
ylabel('$u$')
legend('Numerical method','Analytical','Location','best')
hold off