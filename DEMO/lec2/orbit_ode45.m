% Planetary motion ode45

clear all
close all
clc

global GM;
GM=1;
plot(0,0,'ro','MarkerSize',40, 'MarkerFaceColor', 'r')
tmax=10;  % end time

% Initial condition
x=1; y=0; u=0; v=0.7; 
Y0 = [x; y; u; v;];

hold on
plot(x,y,'bo','MarkerSize',2)

% Use ode45
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
tic
[t, Y] = ode45(@func, [0, tmax], Y0, opts);
toc

plot(Y(:,1),Y(:,2), 'b')
hold off
xlabel('x')
ylabel('y')
title(["Planetary orbit with explicit ode45"])

% Total energy = kinetic energy + potential energy
E = 0.5.*(Y(:,3).^2+Y(:,4).^2) - GM./sqrt(Y(:,1).^2+Y(:,2).^2);

figure(2)
semilogy(t, E)
xlabel('t')
ylabel('E')
title(["Total energy with ode45"])