% Explicit upwind method for 1D linear advection equation with periodic boundary conditions
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
k=2*pi;
u0=sin(k*x);
u=u0; unew=u0; u_ex=u0;

% y-array for computation of exact solution
x=linspace(0,1,jmax);

c=1; % Advection velocity
C=1; % Courant number
% Time step and end of simulation
dt=C*dx/c
tmax=1; % tstop as found in task a)
nmax=ceil(tmax/dt);

% Scale dt stop at tstop
dt=tmax/nmax;
j =[2:jmax];
figure(1)

for n=1:nmax
    
    unew(j) = u(j)*(1-C)+C*u(j-1);
    unew(1) = unew(jmax);
    u=unew;
    t=n*dt;
    u_ex = sin(k*(x-c*t));
    
    plot(x, u, 'x', x, u_ex, 'r')
    %xlim([0 1])
    %ylim([-1,1])
    title(['$t=', num2str(t),'$'])
    
    drawnow;
end

hold on
grid()
xlabel('$y$')
ylabel('$u$')
legend('FTCS','Analytical','Location','best')
hold off

% % Amplification factor of FTCS
% g = 1-4*C*sin(0.5*k*dx)^2;
% u_num = g^nmax*sin(k*x);
% 
% % Amplification factor of exact solution
% g_ex = exp(-k^2*dt);
% u_ex_fourier = g_ex^nmax*sin(k*x); 
% 
% figure(2)
% plot(x, u_num, 'x', x, u_ex_fourier, 'r')
% grid()
% xlabel('$y$')
% ylabel('$u$')
% legend('FTCS from Fourier', 'Exact from Fourier')
% 
% figure(3)
% norm