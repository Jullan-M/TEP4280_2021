% Different methods for 1D linear advection equation with periodic boundary conditions
clear all
close all
clc

set(groot, 'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

% Number of spatial grid points and cell size
jmax=21;
dx=1/(jmax-1);
x=linspace(0,1,jmax)';

% Initialize velocity arrays and set right B.C.
k=4*pi;
u0=cos(k*x);
u=u0; unew=u0; u_ex=u0;

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
j =[2:jmax];
jp =[3:jmax,2]; % Periodic BC u_jmax+1=u_2

% Schemes: 
% 1=FTCS, 2=Expl.Upwind, 3=Impl.Upwind, 4=Lax-Friedrichs, 5=Lax-Wendroff
scheme=5;

if scheme==3
    A=(1+C)*diag(ones(jmax-1,1),0)-C*diag(ones(jmax-2,1), -1);
    A(1,end)=-C;
    b=zeros(jmax-1,1);
end

figure(1)
% Time loop 
for n=1:nmax
    switch scheme
        case 1
            unew(j) =u(j)-0.5*C*(u(jp)-u(j-1));
            % unew(jmax) = u(jmax)-0.5*C*(u(2)-u(jmax-1)); % Periodic BC u_jmax+1=u_2
            unew(1) = unew(jmax);
        case 2
            % Explicit upwind method
            unew(j) = u(j)*(1-C)+C*u(j-1);
            unew(1) = unew(jmax);
        case 3
            % Implicit upwind method
            unew(2:jmax)=A\u(2:jmax);
            unew(1) = unew(jmax);
        case 4
            % Lax-Friedrichs method
            unew(j) = 0.5*(1+C)*u(j-1)+0.5*(1-C)*u(jp);
            unew(1)=unew(jmax);
        case 5
            % Lax-Wendroff method
            unew(j) = u(j)-0.5*C*(u(jp)-u(j-1))+0.5*C^2*(u(jp)-2*u(j)+2*u(j-1));
            unew(1)=unew(jmax);
    end
    u=unew;
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