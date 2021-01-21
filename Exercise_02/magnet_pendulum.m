%% Exercise 2 -- Simulation of chaotic magnetic pendulum

clear all
close all
clc

% Plot magnets
plot(1,0,'ro', -1,0,'ro','MarkerSize',12);

global d; global x1; global x2;
d = 0.1; x1 = 1; x2 = -1;
tmax=10;  % end time
dt=0.001;  % time step size

% Number of time levels
nmax=ceil(tmax/dt); % ceil: rounding up
dt=tmax/nmax;


% Initial condition
t=0;
x=2.0; y=0.1; u=0; v=0; 


hold on
plot(x,y,'bo','MarkerSize',5)

% Loop over time levels with Euler's method
tic
[xp, yp] = euler([x;y;u;v], dt, nmax);
toc
plot(xp,yp, 'b')
hold off
xlabel('x')
ylabel('y')
title(['Chaotic magnetic pendulum, dt=', num2str(dt)])

function [xp, yp] = euler(Y0, dt, nmax)
xp = zeros(1, nmax+1);
yp = xp;
t = 0
Y = Y0;
xp(1)=Y0(1); yp(1)=Y0(1);
for n=1:nmax
    k = func(t, Y);
    Y = Y + dt*k;
    xp(n+1)=Y(1); yp(n+1)=Y(2);
end
end

function [xp, yp] = heun(Y0, dt, nmax)
xp = zeros(1, nmax+1);
yp = xp;
t = 0;
Y = Y0;
xp(1)=Y0(1); yp(1)=Y0(1);
for n=1:nmax
    t = t + dt;
    k1 = func(t, Y);
    k2 = func(t + dt, Y + dt*k1);
    Y = Y + 0.5*(k1+k2)*dt;
    xp(n+1)=Y(1); yp(n+1)=Y(2);
end
end