function dydt = func(t, Y)
global GM
x = Y(1);
y = Y(2);
u = Y(3);
v = Y(4);
r = sqrt(x^2 + y^2);
dydt = [u;
    v;
    -GM*x/r^3;
    -GM*y/r^3];