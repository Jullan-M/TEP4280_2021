% A class that solves ODEs for magnetic pendulum system

classdef ODE_Solver
    properties
        Y0
    end
    methods
        function obj = ODE_Solver(Y0)
            obj.Y0 = Y0'; % Transpose of the initial condition vector
        end
        
        function f = func(obj, t, Y)
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
        end
        
        function [xp, yp, up, vp] = euler(obj, dt, tmax)
            nmax=ceil(tmax/dt);
            xp = zeros(1, nmax+1);
            yp = xp; up = xp; vp = xp;
            t = 0;
            Y = obj.Y0;
            xp(1)=Y(1); yp(1)=Y(2); up(1)=Y(3); vp(1)=Y(4);
            for n=1:nmax
                k = obj.func(t, Y);
                Y = Y + dt*k;
                xp(n+1)=Y(1); yp(n+1)=Y(2); up(n+1)=Y(3); vp(n+1)=Y(4);
            end
        end

        function [xp, yp, up, vp] = heun(obj, dt, tmax)
            nmax=ceil(tmax/dt);
            xp = zeros(1, nmax+1);
            yp = xp; up = xp; vp = xp;
            t = 0;
            Y = obj.Y0;
            xp(1)=obj.Y0(1); yp(1)=obj.Y0(2);
            for n=1:nmax
                t = t + dt;
                k1 = obj.func(t, Y);
                k2 = obj.func(t + dt, Y + dt*k1);
                Y = Y + 0.5*(k1+k2)*dt;
                xp(n+1)=Y(1); yp(n+1)=Y(2); up(n+1)=Y(3); vp(n+1)=Y(4);
            end
        end
    end
end