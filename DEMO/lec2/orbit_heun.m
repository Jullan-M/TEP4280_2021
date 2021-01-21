% Planetary motion Heun's method

clear all
close all
clc

global GM;
GM=1;
plot(0,0,'ro','MarkerSize',40, 'MarkerFaceColor', 'r')
tmax=10;  % end time
dt=0.01;  % time step size

% Number of time levels
nmax=ceil(tmax/dt); % ceil: rounding up
dt=tmax/nmax;

xp = zeros(1, nmax+1);
yp = xp;

t=0;
tim = linspace(0, tmax, nmax+1)
% Initial condition
x=1; y=0; u=0; v=0.7; 
xp(1)=x; yp(1)=y;
Y=[x;y;u;v];
% Total energy = kinetic energy + potential energy
E(1) = 0.5*(u^2+v^2) - GM/sqrt(x^2+y^2);

hold on
plot(x,y,'bo','MarkerSize',2)

% Loop over time levels with Heun's method
tic
for n=1:nmax
    k1 = func(t, Y);
    t = t + dt;
    k2 = func(t + dt, Y + dt*k1);
    Y = Y + 0.5*(k1+k2)*dt;
    xp(n+1)=Y(1); yp(n+1)=Y(2);
    E(n+1) = 0.5*(Y(3)^2+Y(4)^2) - GM/sqrt(Y(1)^2+Y(2)^2);
end
toc
plot(xp,yp, 'b')
hold off
xlabel('x')
ylabel('y')
title(["Planetary orbit with explicit Heun's method, dt=", num2str(dt)])

figure(2)
semilogy(tim, E)
xlabel('t')
ylabel('E')
title(["Total energy with Heun's method, dt=", num2str(dt)])