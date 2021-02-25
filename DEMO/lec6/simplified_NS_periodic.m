% FTCS for 1D heat/diffusion equation with periodic boundary conditions
clear all
close all
clc

set(groot, 'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

% Number of spatial grid points and cell size
jmax=21;
dy=1/(jmax-1);
y=linspace(0,1,jmax);

% Initialize velocity arrays and set right B.C.
k=2*pi;
u0=sin(k*y);
u=u0; unew=u0; u_ex=u0;
U3=zeros(1,jmax);

% y-array for computation of exact solution
y=linspace(0,1,jmax);

% Stable for r <= 0.5
r=0.5 % von Neumann number
% Time step and end of simulation
dt=r*dy^2
tstop=0.05; % tstop as found in task a)
nmax=ceil(tstop/dt);

% Scale dt stop at tstop
dt=tstop/nmax;
j =[2:jmax-1];
figure(1)
for n=1:nmax
    unew(j) = u(j)*(1-2*r)+(u(j-1)+u(j+1))*r;
    unew(jmax) = u(jmax)*(1-2*r)+(u(jmax-1)+u(2))*r;
    u=unew;
    t=n*dt;
    u_ex = exp(-k^2*t)*u0;
    
    plot(y, u, 'x', y, u_ex, 'r')
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

% Amplification factor of FTCS
g = 1-4*r*sin(0.5*k*dy)^2;
u_num = g^nmax*sin(k*y);

% Amplification factor of exact solution
g_ex = exp(-k^2*dt);
u_ex_fourier = g_ex^nmax*sin(k*y); 

figure(2)
plot(y, u_num, 'x', y, u_ex_fourier, 'r')
grid()
xlabel('$y$')
ylabel('$u$')
legend('FTCS from Fourier', 'Exact from Fourier')

figure(3)
norm