% Stiff ODE example: Flame propagation model ODE dy/dt = y^2 - y^3
clear all
close all
clc

delta = 0.01;
tend = 2/delta
dt =  0.01;
nmax = ceil(tend/dt);
y = delta;

% RHS of ODE defined as an anonymous function
f = @(t,y) y^2 - y^3;

opts = odeset('Stats', 'on', 'RelTol',1.e-4)
tic
[t, y] = ode23s(f, [0, tend], delta, opts);
toc

plot(t,y, 'x')
hold on
a = 1/delta -1;
ya = 1./(lambertw(a*exp(a-t)) + 1)
plot(t, ya, 'r')
xlabel('t')
ylabel('y')
title('Flame propagation model ODE dy/dt = y^2 - y^3 with \delta = ', num2str(delta))
