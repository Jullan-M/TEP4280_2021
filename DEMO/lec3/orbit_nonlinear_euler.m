% Planetary motion with implicit Euler's method
% solving the nonlinear system of equations each time step

clear all
close all
clc

global GM;
GM=1;

tmax=10;  % end time
dt=0.001;  % time step size

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

opts = optimoptions('fsolve', 'Display', 'none');

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
    % define nonlinear function Y -dt*f(t_n+1, Y) - Y_n = 0 where Y=Y_n+1
    xi = Yi(1); yi = Yi(2); ui = Yi(3); vi = Yi(4);
    res = @(z) [z(1)-dt*z(3)-xi;
                z(2)-dt*z(4)-yi;
                z(3)+dt*GM*z(1)/sqrt(z(1)^2+z(2)^2)^3 - ui;
                z(4)+dt*GM*z(2)/sqrt(z(1)^2+z(2)^2)^3 - vi];
    Yi = fsolve(res, Yi, opts);
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