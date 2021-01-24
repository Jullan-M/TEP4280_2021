clear all
close all
clc

global d; global x1; global x2;
d = 0.1; x1 = 1; x2 = -1;

sol1 = ODE_Solver([2.0,0.1000,0,0]);
sol2 = ODE_Solver([2.0,0.10001,0,0]);

tmax=10;  % end time

% Solve ODE with ODE45
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
tic
[t1, Y1] = ode45(@sol1.func, [0, tmax], [2.0,0.1000,0,0]);
toc
tic
[t2, Y2] = ode45(@sol2.func, [0, tmax], [2.0,0.10001,0,0]);
toc

figure(1)
hold on
plot(Y1(:,1), Y1(:,2), 'r', Y2(:,1), Y2(:,2), 'b') % Plot results
plot(2.0,0.1000,'o','MarkerSize',5); % Plot start position
plot(1,0,'ro', -1,0,'ro','MarkerSize',12); % Plot magnets
hold off
legend(["$y=0.10000$", "$y=0.10001$"], 'Location', 'best')
xlabel('x')
ylabel('y')
title(['Chaotic magnetic pendulum, ODE45'])