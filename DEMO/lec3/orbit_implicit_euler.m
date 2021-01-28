% Planetary motion with implicit Euler's method
% approximating r^3_n+1 by r^3_n

clear all
close all
clc

global GM;
GM=1;

tmax=10;  % end time
dt=0.00001;  % time step size

% Number of time levels
nmax=ceil(tmax/dt); % ceil: rounding up
dt=tmax/nmax;

xp = zeros(1, nmax+1);
yp = xp;

t=0;
% Initial condition
x=1; y=0; u=0; v=0.7; xp(1)=x; yp(1)=y;
Yi = [x;y;u;v];
YY(:,1) = Yi;

% Define A matrix. -1 values will be replaced by dt*GM/r3 in the loop
A = [1 0 -dt 0;
    0 1 0 -dt;
    -1 0 1 0;
    0 -1 0 1];

% Loop over time levels with both the explicit and implicit Euler's method
tic
for n=1:nmax
    % Explicit Euler
    xn=x+dt*u;
    yn=y+dt*v;
    r3=(sqrt(x^2+y^2))^3;
    u=u-dt*GM*x/r3;
    v=v-dt*GM*y/r3;
    x=xn;
    y=yn;
    xp(n+1)=x; yp(n+1)=y;
    
    % Implicit Euler
    xi = Yi(1); yi = Yi(2);
    r3=(sqrt(xi^2+yi^2))^3;
    gv = dt*GM/r3;
    A(3,1) = gv; A(4,2) = gv;
    Yi = A \ Yi;
    YY(:,n+1) = Yi;
end
toc

hold on
plot(xp,yp, 'b', YY(1,:), YY(2,:), 'g')
plot(0,0,'ro','MarkerSize',40)
plot(xp(1),yp(1),'mo','MarkerSize',3)
hold off
xlabel('x')
ylabel('y')
legend({'Explicit Euler method', 'Implicit Euler method'}, 'Location', 'best')
title(['Planetary orbit with dt=', num2str(dt)])