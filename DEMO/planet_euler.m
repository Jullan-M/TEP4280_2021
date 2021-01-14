clear all
close all
clc
% 1st order forward Euler method

GM = 1;
dt = 0.01;
tmax = 100;
x0 = 1; y0 = 0; u0 = 0; v0 = 1;
x = x0; y = y0; u = u0; v = v0;

figure(1);
h=plot(x, y, 'bo');
axis([-1 1 -1 1]);

hold on
plot(0,0, 'ro');
for i=0:round(tmax/dt)
    r = sqrt(x^2 + y^2);
    xn = x + dt*u;
    yn = y + dt*v;
    u = u + dt*(-GM*x/r^3);
    v = v + dt*(-GM*y/r^3);
    x=xn; y=yn;
    %plot(x,y, 'bo')
    set(h, 'XData', x, 'YData', y);
    if mod(i,100) == 0
        plot(x,y,'b.')
    end
    drawnow;
end
hold off
