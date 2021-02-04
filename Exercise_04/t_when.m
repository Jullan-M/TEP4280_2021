clear all
close all
clc

un = @(t, y, i) erfc((2*i+1-y)/(2*sqrt(t)))-erfc((2*i+1+y)/(2*sqrt(t)));
u = @(t, y, i) sum(un(t, y, i));

i = (0:5);
tstop = fzero(@(t) u(t, 0.5, i)-0.25, [0 1]);
num2str(tstop, 30)