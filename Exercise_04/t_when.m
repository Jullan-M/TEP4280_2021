% Task a) Finding zeros of f(y,t) = u(y=0.5, t) - 0.25 = 0
clear all
close all
clc

% Element i of the sum in equation (3).
un = @(t, y, i) erfc((2*i+1-y)/(2*sqrt(t)))-erfc((2*i+1+y)/(2*sqrt(t)));
% Equation (3) summed up to index n
u = @(t, y, n) sum(un(t, y, 0:n));

n = 5;
tstop = fzero(@(t) u(t, 0.5, n)-0.25, [0 1]);
num2str(tstop, 30)