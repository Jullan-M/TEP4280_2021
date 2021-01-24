clear all
close all
clc

global d; global x1; global x2;
d = 0.1; x1 = 1; x2 = -1;
x=2.0; y=0.1000; u=0; v=0;

solver = ODE_Solver([x,y,u,v]);

tmax=10;  % end time
dts = 0.000025:0.0000025:0.00005;

x_ex = 0.6275405986573285854746018230798654258251190185546875;
p_h = zeros(6,1); p_e = zeros(6,1);
for i=1:11
    dt = dts(i)
    [x2t_h, y2t_h] = solver.heun(2*dt, tmax);
    [xt_h, yt_h] = solver.heun(dt, tmax);
    p_h(i) = log(abs((x2t_h(end)- x_ex)/(xt_h(end) - x_ex)))/log(2);
    [x2t_e, y2t_e] = solver.euler(2*dt, tmax);
    [xt_e, yt_e] = solver.euler(dt, tmax);
    p_e(i) = log(abs((x2t_e(end)- x_ex)/(xt_e(end) - x_ex)))/log(2);
end


figure(4)
hold on
plot(dts, p_e, 'b', dts, p_h, 'r') % Plot results
hold off
legend(["Euler's method", "Heun's method"], 'Location', 'best')
xlabel('dt')
ylabel('p')
title(['Chaoticmagnetic pendulum'])