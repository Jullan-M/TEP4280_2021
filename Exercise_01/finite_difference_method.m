clear all
close all
clc

tan1_ex = tan(1)^2 + 1;
n = 1000
h = logspace(-16, -1, n);

%% 1st order forward difference with double precision
tan_I = ((tan(1+h) - tan(1)))./h;


%% 2nd order central difference with double precision
tan_II = ((tan(1+h) - tan(1-h)))./(2*h);

%% 2nd order finite difference with single precision
tan_III = ((single(tan(1+h)) - single(tan(1)) ))./(2*single(h) );

% Plot curves
loglog(1,1,1,1,1,1,h,tan_I-tan1_ex,'--',h,tan_II-tan1_ex,'--',h,tan_III-tan1_ex,'--')
legend({'1. order forward or backward difference double precision ','2. order central difference , double precision','2. order central difference , single precision','1. order','2. order' ,'-1. order'},'Location','best')
% Use a TeX-interpreter to get a nicer title
title('$\mathrm{d}(\tan x)/ \mathrm{d}x = \tan^2(x)+1$','Interpreter','latex')
xlabel('h')
ylabel('err= |(df/dx_{FD} - df/dx_{An})/df/dx_{An}|')
grid on
