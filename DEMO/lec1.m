%% Machine epsilon
clear all % Clears all data in memory
close all % Closes all open figure windows
clc % Clear command window
eps % Built-in value, smallest value that can make a difference to 1

%% Double precision
eps1 = 1;
while (eps1+1) > 1
    eps1 = eps1/2;
end
eps1 = 2*eps1

%% Single precision
eps2 = single(1);
while (eps2+1) > 1
    eps2 = eps2/2;
end
eps2 = 2*eps2

plot(1,1)