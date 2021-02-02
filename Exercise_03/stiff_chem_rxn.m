clear all
close all
clc

k1=0.04; k2=3e7; k3=1e4;

y0 = [1;0;0;];
tend = 0.3;

% Task 2 a)
f = @(t, y) [-k1*y(1)+k3*y(2)*y(3);
    k1*y(1)-k2*y(2)^2-k3*y(2)*y(3);
    k2*y(2)^2];
% Solve ODEs using different in-built MATLAB methods
[t_45, Y_45] = ode45(f, [0, tend], y0);
[t_15s, Y_15s] = ode15s(f, [0, tend], y0);
[t_23s, Y_23s] = ode23s(f, [0, tend], y0);

% Task 2 b)
set(groot, 'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
% Plot y2 computed with ode45, ode15s and ode23s in the same figure
figure(1)
hold on
plot(t_45,Y_45(:,2), 'r', t_15s,Y_15s(:,2), 'g', t_23s,Y_23s(:,2), 'b--')
hold off
legend(["ODE45", "ODE15S", "ODE23S"])
ylabel("$y_2(t)$")
xlabel("$t$")
grid()

% Plot |y1+y2+y3-1|
figure(2)
semilogy(t_45, abs(Y_45(:,1)+Y_45(:,2)+Y_45(:,3)- 1), 'r')
hold on
semilogy(t_15s, abs(Y_15s(:,1)+Y_15s(:,2)+Y_15s(:,3)- 1), 'g')
semilogy(t_23s, abs(Y_23s(:,1)+Y_23s(:,2)+Y_23s(:,3)- 1), 'b')
hold off
legend(["ODE45", "ODE15S", "ODE23S"])
ylabel("$|y_1+y_2+y_3-1|$")
xlabel("$t$")
ylabel("$|y_1+y_2+y_3-1|(t)$")
grid()

% Task 2 c)
tend = 1e11;
% Set relative and absolute tolerance lower for higher accuracy
opts = odeset('RelTol' ,1e-7, 'AbsTol',1e-12);
[t, Y] = ode23s(f, [0, tend], y0, opts);

figure(3)
semilogx(t, Y(:,1), 'r')
ylabel("$y_1(t)$")
xlabel("$t$")
title("Species A")
grid()

figure(4)
semilogx(t, Y(:,2), 'g')
ylabel("$y_2(t)$")
xlabel("$t$")
title("Species B")
grid()

figure(5)
semilogx(t, Y(:,3), 'b')
ylabel("$y_3(t)$")
xlabel("$t$")
title("Species C")
grid()

% Plot sum of species in the same figure (should be equal to 1)
figure(6)
semilogx(t, Y(:,1)+Y(:,2)+Y(:,3), 'b')
ylabel("$y_1(t)+y_2(t)+y_3(t)$")
xlabel("$t$")
title("Sum of Species A, B and C")
grid()