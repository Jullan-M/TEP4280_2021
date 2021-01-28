clear all
close all
clc

global d; global x1; global x2;
d = 0.1; x1 = 1; x2 = -1;
x=2.0; y=0.1000; u=0; v=0;

solver = ODE_Solver([x,y,u,v]);

tmax=10;  % end time
dts = 0.0000125:0.00000625:0.00005;
n = length(dts);

% "Exact" value of the x-coordinate of the end position (Heun's method with dt=0.000001)
xt_ex = 0.62754030005734617869705971315852366387844085693359;
p_h = zeros(n,1); p_e = zeros(n,1);
for i=1:n
    dt = dts(i)
    [x2t, y2t] = solver.heun(2*dt, tmax);
    [xt, yt] = solver.heun(dt, tmax);
    % Calculate the end position error and the convergence rate
    p_h(i) = log(abs((x2t(end) - xt_ex)/(xt(end) - xt_ex)))/log(2);
    
    [x2t, y2t] = solver.euler(2*dt, tmax);
    [xt, yt] = solver.euler(dt, tmax);
    p_e(i) = log(abs((x2t(end) - xt_ex)/(xt(end) - xt_ex)))/log(2);
end

set(groot, 'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');
figure(1)
hold on
plot(2*dts, p_e, 'b', 2*dts, p_h, 'r') % Plot results
hold off
legend(["Euler's method", "Heun's method"], 'Location', 'best')
xlabel('$\Delta t$')
ylabel('$p$')
title(['Convergence rate of numerical methods'])
grid()