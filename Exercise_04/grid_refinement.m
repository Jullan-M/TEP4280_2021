% Grid refinement study of FTCS scheme for the velocity profile problem
clear all
close all
clc

mmax=8; % max number of iterations of grid refinement
rs=[0.5, 1/6];
tstop=0.0946869595678490744150579416782; % tstop as found in task a)

relerror = zeros(2,mmax);
gridpts = 2.^(1:mmax)+1;

% Exact solution at u(y=0.5, t=tstop) (see task a)
u_ex= 0.25; % Analytic solution

for l=1:2
    r = rs(l);
    for m=1:mmax
        m
        % Number of spatial grid points and cell size
        jmax=2^m+1;
        dy=1/(jmax-1);
        mid=ceil(jmax/2); % Middle point y=1/2

        % Initialize velocity arrays and set right B.C.
        u=zeros(jmax,1);
        U3=zeros(1,jmax);
        u(jmax)=1;

        % y-array for computation of exact solution
        y=linspace(0,1,jmax);    

        % Adjust time step such that r=0.5
        dt=r*dy^2;
        nmax=ceil(tstop/dt);

        for n=1:nmax
            u(2:jmax-1) = u(2:jmax-1)*(1-2*r)+(u(3:jmax)+u(1:jmax-2))*r;
            t=dt*n;
        end
        relerror(l,m) = abs((u(mid)-u_ex)/u_ex);
    end
end
set(groot, 'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

% Plot the relative error
loglog(gridpts,relerror(1,:),'r-.', gridpts, relerror(2,:),'b-.')
legend('$r=0.5$','$r=1/6$')
hold on
xlabel('Grid points, $j_{max}$')
ylabel('Relative error, $|(u(0.5, t_{s})-u_{exact}(0.5, t_{s}))/u_{exact}(0.5, t_{s})|$')

grid()
hold off