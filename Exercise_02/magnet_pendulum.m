%% Exercise 2 -- Simulation of chaotic magnetic pendulum
clear all
close all
clc

global d; global x1; global x2;
d = 0.1; x1 = 1; x2 = -1;

% Initial conditions
x=2.0; y=0.1000; u=0; v=0;

solver = ODE_Solver([x,y,u,v]);

tmax=10;  % end time

% Solve ODE with Euler's method, Heun's method and ODE45 with dt = 0.001
dt1=0.001;  % time step size
tic
[xp_eul, yp_eul, up_eul, vp_eul] = solver.euler(dt1, tmax);
toc

tic
[xp_heu, yp_heu, up_heu, vp_heu] = solver.heun(dt1, tmax);
toc

tic
opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
[t_45, Y_45] = ode45(@solver.func, [0, tmax], [x,y,u,v], opts);
toc

set(groot, 'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
figure(1)
hold on
plot(xp_eul,yp_eul, 'b', xp_heu,yp_heu, 'g', Y_45(:,1), Y_45(:,2), '--m') % Plot results
plot(x,y,'o','MarkerSize',5); % Plot start position
plot(1,0,'ro', -1,0,'ro','MarkerSize',12); % Plot magnets
hold off
legend(["Euler's method", "Heun's method", "ODE45"], 'Location', 'best')
xlabel('$x$')
ylabel('$y$')
grid()
title(['Chaotic magnetic pendulum, $\Delta t=', num2str(dt1), '$'])

ts = 0:0.001:tmax;
% Calculate energies for each of the numerical methods
E_eul = 0.5*(up_eul.^2 + vp_eul.^2) - (1./sqrt((xp_eul-x1).^2+yp_eul.^2 + d.^2) + 1./sqrt((xp_eul-x2).^2+yp_eul.^2 + d.^2));
E_heu = 0.5*(up_heu.^2 + vp_heu.^2) - (1./sqrt((xp_heu-x1).^2+yp_heu.^2 + d.^2) + 1./sqrt((xp_heu-x2).^2+yp_heu.^2 + d.^2));
E_o45 = 0.5*(Y_45(:,3).^2 + Y_45(:,4).^2) - (1./sqrt((Y_45(:,1)-x1).^2+Y_45(:,2).^2 + d.^2) + 1./sqrt((Y_45(:,1)-x2).^2+Y_45(:,2).^2 + d.^2));

% Plot mechanical energy with dt = 0.001
figure(4)
hold on
subplot(3,1,1)
plot(ts,E_eul, 'b')
title(["$E_{mech}$ for Euler's method"])
xlabel('$t$')
ylabel('$E$')
grid()

subplot(3,1,2)
plot(ts,E_heu, 'g')
title(["$E_{mech}$ for Heun's method"])
xlabel('$t$')
ylabel('$E$')
grid()

subplot(3,1,3)
plot(t_45,E_o45, 'm')
title(["$E_{mech}$ for ODE45"])
hold off
xlabel('$t$')
ylabel('$E$')
grid()

% Solve ODE with Euler's method, Heun's method with dt = 0.0001
dt2=0.0001;
tic
[xp_eul, yp_eul] = solver.euler(dt2, tmax);
toc

tic
[xp_heu, yp_heu] = solver.heun(dt2, tmax);
toc

figure(2)
hold on
plot(xp_eul,yp_eul, 'b', xp_heu,yp_heu, 'g', Y_45(:,1), Y_45(:,2), '--m') % Plot results
plot(x,y,'o','MarkerSize',5); % Plot start position
plot(1,0,'ro', -1,0,'ro','MarkerSize',12); % Plot magnets
hold off
legend(["Euler's method", "Heun's method", "ODE45"], 'Location', 'best')
xlabel('$x$')
ylabel('$y$')
grid()
title(['Chaotic magnetic pendulum, $\Delta t=', num2str(dt2), '$'])


% Solve ODE with Euler's method, Heun's method with dt = 0.0001
dt3=0.00001;
tic
[xp_eul, yp_eul] = solver.euler(dt3, tmax);
toc

tic
[xp_heu, yp_heu] = solver.heun(dt3, tmax);
toc

figure(3)
hold on
plot(xp_eul,yp_eul, 'b', xp_heu,yp_heu, 'g', Y_45(:,1), Y_45(:,2), '--m') % Plot results
plot(x,y,'o','MarkerSize',5); % Plot start position
plot(1,0,'ro', -1,0,'ro','MarkerSize',12); % Plot magnets
hold off
legend(["Euler's method", "Heun's method", "ODE45"], 'Location', 'best')
xlabel('$x$')
ylabel('$y$')
grid()
title(['Chaotic magnetic pendulum, $\Delta t=', num2str(dt3), '$'])