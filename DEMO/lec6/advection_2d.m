% Explicit upwind method for 2D linear advection equation 
% u_t + c_1 u_x + c_2 u_y = 0, where c_1 and c_2 are the
% advection velocity components in x- and y-directions, resp.
% Here, c_1 = c_2 and dx = dy are assumed, i.e., 
% transport in north-east direction on a Cartesian grid.

clear all
close all
clc

a=10;
b=7;

imax=101; %401; 201; 101;
jmax= 71; %281; 141;  71;
xa=linspace(-a,a,imax);
h=xa(2)-xa(1);           % h=dx=xy
ya=-b+[0:jmax-1]*h;
[x,y]=meshgrid(xa,ya);
x=x'; y=y';

u=x*0;

C=0.5; % Courant number C=c*dt/h
% Stability bound: 0 <= C=c*dt/h <= 1/2 for h=dx=dy and c=c_1=c_2
% General stability bound: 
% c_1*dt/dx + c_2*dt/dy <= 1, 0 <= c_1*dt/dx, 0 <= c_2*dt/dy

c=1; % c=c_1=c_2
dt=C*h/c;
tmax=5; 
nmax=ceil(tmax/dt);
dt=tmax/nmax;
i=[2:imax];
j=[2:jmax];
C=c*dt/h;
A=1;

u=A*exp(-(x.^2+y.^2));   % Gauss pulse at origin

nu=u;

surf(x,y,u)
pause
t=0;

for n=1:nmax
    % explicit upwind FDM
    nu(i,j)=u(i,j)-C*(u(i,j)-u(i-1,j)+u(i,j)-u(i,j-1));
    % nu(i,j)=u(i-1,j-1); % gives exact solution for C=1
    
    u=nu;
    t=t+dt;
    surf(x,y,u)
    colormap('jet');
    title(['t=',num2str(t,2)])
    zlim([0 1])
    drawnow
end

xlabel('x')
ylabel('y')
zlabel('u')
title(['Explicit upwind scheme for 2D linear advection equation,t=',num2str(t),', C=',num2str(C),', imax=',num2str(imax),', jmax=',num2str(jmax)])
colorbar

figure(2)
u_exact=A*exp(-((x-c*t).^2+(y-c*t).^2));
surf(x,y,u_exact)
colormap('jet');
xlabel('x')
ylabel('y')
zlabel('u_{exact}')
title(['Exact soluition of 2D linear advection equation,t=',num2str(t),', imax=',num2str(imax),', jmax=',num2str(jmax)])
zlim([0 1])
colorbar

disp(['||u - u_exact||_2 = ',num2str(h*norm(u-u_exact))])