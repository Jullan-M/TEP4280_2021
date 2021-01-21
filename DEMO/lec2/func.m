% RHS of ODE for planetary orbit around the sun

function f = func(t, Y)
global GM;
f=zeros(4,1);
x = Y(1); y = Y(2); u = Y(3); v = Y(4);

f(1) = u; f(2) = v;

r3 = (sqrt(x^2+y^2))^3;
f(3) = -GM*x/r3;
f(4) = -GM*y/r3;