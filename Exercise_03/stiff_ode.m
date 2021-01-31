clear all
close all
clc

dt = 0.01;
y0 = 10;
tend = 10;

sol = ODE_Solver(y0, @func);
y_eeul = sol.exp_euler(dt, tend);
y_ieul = sol.imp_euler(dt, tend);
t = 0:dt:tend;
y_exact = 1 + t + (y0 - 1)*exp(-100.*t);
plot(t, y_exact, 'b', t, y_eeul, 'r--', t, y_ieul, 'm--')
legend(["Exact solution", "Explicit Euler", "Implicit Euler"])
xlabel("t")
ylabel("y(t)")

function dydt = func(t, y)
    dydt = -100*y + 100*t + 101;
end