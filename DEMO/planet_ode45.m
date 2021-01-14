clear all
close all
clc
% 1st order forward Euler method
global GM
GM = 1;
dt = 0.01;
tmax = 100;
x0 = 1; y0 = 0; u0 = 0; v0 = 1;
x = x0; y = y0; u = u0; v = v0;

figure(1);
h=plot(x, y, 'bo');
axis([-1 1 -1 1]);

options = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t, y] = ode45(@func, [0 tmax], [x y u v], options);
plot(0, 0, 'ro')
plot(y(:,1), y(:,2))