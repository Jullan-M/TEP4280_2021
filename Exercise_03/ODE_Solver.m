% A class that solves ODEs for magnetic pendulum system given a set of ICs

classdef ODE_Solver
    properties
        Y0
        func
    end
    methods
        function obj = ODE_Solver(Y0, func)
            obj.Y0 = Y0;
            obj.func = func;
        end
        
        % Explicit Euler method
        function y = exp_euler(obj, dt, tmax)
            nmax=ceil(tmax/dt);
            y = zeros(1, nmax+1);
            t = 0;
            Y = obj.Y0;
            y(1) = Y;
            for n=1:nmax
                t = t + dt;
                k = obj.func(t, Y);
                Y = Y + dt*k;
                y(n+1) = Y;
            end
        end
        
        % Implicit Euler method
        function y = imp_euler(obj, dt, tmax)
            nmax=ceil(tmax/dt);
            y = zeros(1, nmax+1);
            dfdy = obj.func(0, 1) - obj.func(0, 0);
            t = 0;
            Y = obj.Y0;
            y(1) = Y;
            for n=1:nmax
                t = t + dt;
                % Only one iteration of Newton's method is applied here
                dy = obj.func(t, Y) / (1 - dt*dfdy);
                Y = Y + dt*dy;
                y(n+1) = Y;
            end
        end
        
        % TODO:
        % Trapezoidal method
        function y = trapezoid(obj, dt, tmax)
            nmax=ceil(tmax/dt);
            y = zeros(1, nmax+1);
            dfdy = obj.func(0, 1) - obj.func(0, 0);
            t = 0;
            Y = obj.Y0;
            y(1) = Y;
            for n=1:nmax
                t = t + dt;
                % Only one iteration of Newton's method is applied here
                dy = obj.func(t, Y) / (1 - 0.5*dt*dfdy);
                Y = Y + dt*dy;
                y(n+1) = Y;
            end
        end
    end
end