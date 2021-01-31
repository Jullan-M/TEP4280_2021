clear all
close all
clc

dt = 0.01;
y0 = 10;
tend = 10;

sol = ODE_Solver(y0, @func);
y_eeul = sol.exp_euler(dt, tend);
y_ieul = sol.imp_euler(dt, tend);
y_trap = sol.trapezoid(dt, tend);

t = 0:dt:tend;
y_exact = 1 + t + (y0 - 1)*exp(-100.*t);

figure(1)
plot(t, y_exact, 'k', t, y_eeul, 'r--', t, y_ieul, 'm--', t, y_trap, 'b--')
legend(["Exact solution", "Explicit Euler", "Implicit Euler", "Trapezoidal"])
xlabel("t")
ylabel("y(t)")
grid()

figure(2)
title(["Error"])
plot(t, abs(y_eeul-y_exact), 'r--', t, abs(y_ieul-y_exact), 'm--', t, abs(y_trap-y_exact), 'b--')
legend(["Explicit Euler", "Implicit Euler", "Trapezoidal"])
xlabel("t")
ylabel("y(t)")
grid()

function dydt = func(t, y)
    dydt = -100*y + 100*t + 101;
end