%% Exercise 1 -- TEP4280
%%
clear all
close all
clc

%% a) Calculating the derivative using the finite difference method
tan1_ex = tan(1)^2 + 1; % Exact value of dtan(x)/dx  at x=1
n = 1000;
h1 = logspace(-16, -1, n);
h2 = logspace(-9, -1, n);

% 1st order forward difference with double precision
tan_I = ((tan(1+h1) - tan(1)))./h1;

% 2nd order central difference with double precision
tan_II = ((tan(1+h1) - tan(1-h1)))./(2*h1);

% 2nd order central difference with single precision
tan_III = ((tan(1+single(h1)) - tan(1-single(h1)) ))./(2*single(h1) );

% 3rd order forward difference with double precision
tan_IV = (-11/6*tan(1) + 3*tan(1+h2) -3/2*tan(1+2*h2) +1/3*tan(1+3*h2) )./h2;

% 4th order central difference with double precision
tan_V = (1/12*tan(1-2*h2) - 2/3*tan(1-h2) + 2/3*tan(1+h2)-1/12*tan(1+2*h2) )./(h2);

% 4th order central difference with single precision
tan_VI = (1/12*tan(1-2*single(h2)) - 2/3*tan(1-single(h2)) + 2/3*tan(1+single(h2))-1/12*tan(1+2*single(h2)) )./(single(h2));


err_I = abs(tan_I-tan1_ex)/tan1_ex;
err_II = abs(tan_II-tan1_ex)/tan1_ex;
err_III = abs(tan_III-tan1_ex)/tan1_ex;
err_IV = abs(tan_IV-tan1_ex)/tan1_ex;
err_V = abs(tan_V-tan1_ex)/tan1_ex;
err_VI = abs(tan_VI-tan1_ex)/tan1_ex;

%% b) First and second order FDM
% Here we see that while the trunctuation error scales in accordance with
% the order of the FDM method used as the the grid spacing $\mathrm{d}x$
% gets smaller, it will inevatibly reach a point where rounding errors
% become very significant.
% 
% In fact, for a grid spacing at around the order of $10^{-8}$, the 2nd order
% method for double precision numbers is no better than the 1st order
% method for double precision numbers. The same holds true for the single
% precision numbers, except that their point of dimishing returns is
% even larger at $\mathrm{d}x \sim 10^{-3}$. It is also noteworthy that the
% error becomes even larger for $\mathrm{d}x$ smaller than this point.

% Use a TeX-interpreter to get nicer text
set(groot, 'defaultTextInterpreter','latex');

figure(1)
loglog(h1,err_I,'--',h1,err_II,'--',h1,err_III,'--')
legend({'1st order forward difference, double precision ','2nd order central difference, double precision','2nd order central difference, single precision'},'Location','best')

title('$\mathrm{d}(\tan x)/ \mathrm{d}x = \tan^2(x)+1$')
xlabel('$h$')
ylabel('error$= \left|(\mathrm{d}f/\mathrm{d}x_{FD} - \mathrm{d}f/\mathrm{d}x_{An})/\mathrm{d}f/\mathrm{d}x_{An}\right|$')
grid on

%% c) Third and fourth order FDM
% The same principle applies for higher order finite difference methods.
% For double presision numbers, the point where rounding errors become
% apparent is at around $10^{-5}$ and $10^{-4}$ for 3rd and 4th order FDM,
% respectively. For single precision numbers, this point is as high as
% $\sim 10^{-2}$.
%
% The takeaway is that while there is precision to be gained in decreasing
% the grid spacing in FDM, there will come a point where one should be
% mindful of rounding errors that will nullify this precision. Where this
% point is depends both on the precision of the number format used and on
% the order of finite difference method.
figure(2)
loglog(h2,err_IV,'--', h2,err_V,'--', h2,err_VI,'--', h2, h2.^3, '-', h2, h2.^4, '-')
ylim([10^(-15), inf])
legend({'3rd order fwd diff., DP', '4th order cent diff., DP', '4th order cent diff., SP', '3rd order slope', '4th order slope'}, 'Location', 'best')
title('$\mathrm{d}(\tan x)/ \mathrm{d}x = \tan^2(x)+1$')
xlabel('$h$')
ylabel('error$= \left|(\mathrm{d}f/\mathrm{d}x_{FD} - \mathrm{d}f/\mathrm{d}x_{An})/\mathrm{d}f/\mathrm{d}x_{An}\right|$')
grid on
