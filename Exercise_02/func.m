% RHS of ODEs for magnetic pendulum system

function f = func(t, Y)
global d; global x1; global x2;
f=zeros(4,1);
x = Y(1); y = Y(2); u = Y(3); v = Y(4);

yz2 = y^2 + d^2;
r3_I=(sqrt((x-x1)^2+yz2))^3;
r3_II=(sqrt((x-x2)^2+yz2))^3;

f(1) = u; 
f(2) = v;
f(3) = -((x-x1)/r3_I + (x-x2)/r3_II);
f(4) = -y*(1/r3_I + 1/r3_II);