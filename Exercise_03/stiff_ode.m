clear all
close all
clc

dts = [0.01, 0.1, 0.25];
y0 = 10;
tend = 10;

% Define differential function and ODE_Solver object
func = @(t,y) -100*y + 100*t + 101;
sol = ODE_Solver(y0, func);
t_ex = 0:0.01:tend;
y_exact = 1 + t_ex + (y0 - 1)*exp(-100.*t_ex);

set(groot, 'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaultLegendInterpreter','latex');

for i=1:3
    %y_eeul = sol.exp_euler(dt, tend);
    %y_ieul = sol.imp_euler(dt, tend);
    dt = dts(i);
    y_trap = sol.trapezoid(dt, tend);
    t = 0:dt:tend;
    
    
    figure(i)
    plot(t_ex, y_exact, 'k', t, y_trap, 'r--')
    legend(["Exact solution", "Trapezoidal method"])
    xlabel("$t$")
    ylabel("$y(t)$")
    title(['$\Delta t =', num2str(dt), '$']);
    grid()
    
    
    error = abs(y_trap-y_exact(1:dt/0.01:end));
    max(error)
end