clear all
close all
clc

kdx = linspace(0,pi,201); % kdx from 0 to pi

% Amplification factor as a function of k dx and r
g_ftcs = @(kdx, r) 1-4*r*sin(kdx./2).^2; % FTCS scheme
g_btcs = @(kdx, r) 1./(1+4*r*sin(kdx./2).^2); % Simple implicit method
g_ex = @(kdx, r) exp(-r*(kdx).^2); % Exact solution

set(groot, 'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

r=0.5; % Plotting r = 0.5
figure(1)
plot(kdx, g_ftcs(kdx, r), kdx, g_btcs(kdx, r), kdx, g_ex(kdx, r), '--') 
title(['von Neumann number, $r=', num2str(r), '$'])
legend(["FTCS scheme", "BTCS scheme", "Exact solution"])
ylabel("Amplification factor, $g(k \Delta x)$")
xlabel("$k \Delta x$")
grid()

r=1/6; % Plotting r = 1/6
figure(2)
plot(kdx, g_ftcs(kdx, r), kdx, g_btcs(kdx, r), kdx, g_ex(kdx, r), '--')
title(['von Neumann number, $r=', num2str(r), '$'])
legend(["FTCS scheme", "BTCS scheme", "Exact solution"])
ylabel("Amplification factor, $g(k \Delta x)$")
xlabel("$k \Delta x$")
grid()