% Planetary motion with explicit Euler's method

clear all
close all
clc

GM=1;
plot(0,0,'ro','MarkerSize',40)
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

hold on
plot(x,y,'bo','MarkerSize',2)

% Loop over time levels with Euler's method
tic
for n=1:nmax
    xn=x+dt*u;
    yn=y+dt*v;
    r3=(sqrt(x^2+y^2))^3;
    u=u-dt*GM*x/r3;
    v=v-dt*GM*y/r3;
    x=xn;
    y=yn;
    xp(n+1)=x; yp(n+1)=y;
end
toc
plot(xp,yp, 'b')
hold off
xlabel('x')
ylabel('y')
title(['Planetary orbit with explicit Euler method, dt=', num2str(dt)])