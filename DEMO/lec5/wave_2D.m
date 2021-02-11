% CTCS scheme for 2D wave equation u_tt = c^2 ( u_xx + u_yy),
% where c is the speed of sound.
% The elliptic domain is supposed to be open to the ambient such the
% the acoustic pressure outside the elliptic domain is zero.


clear all
close all
clc

a=10;
b=7;

imax=201;
jmax=141;
xa=linspace(-a,a,imax);
h=xa(2)-xa(1);           % h=dx=dy
ya=-b+[0:jmax-1]*h; 
[x,y]=meshgrid(xa,ya);
x=x'; y=y';

u=x*0;

c=1;
dt=h/(c*sqrt(2)); % Courant number C = c*dt/h = 1/sqrt(2) at stability limit
tmax=50;
nmax=ceil(tmax/dt);
dt=tmax/nmax;
i=[2:imax-1];
j=[2:jmax-1];
C=c*dt/h;         % Courant number C = c*dt/h calculated with current dt.
A=1;

% u(28:34,68:74)=A;                      % square pulse at one focal point
u=A*exp(-((x+sqrt(a^2-b^2)).^2+y.^2)); % Gauss pulse at one focal point

u(imax,jmax)=-A;          % for scaling the plot
u(1,jmax)=A;              % for scaling the plot
nu=u;
ou=u;

surf(x,y,u)
pause
t=0;
e=(x/a).^2+(y/b).^2<1; % e(x,y)=1, if (x,y) in ellipse; e(x,y)=0 else

for n=1:nmax
    % CTCS FDM for dx=dy
    nu(i,j)=2*u(i,j)-ou(i,j)+C^2*((u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1))-4*u(i,j));
    nu(i,j)=nu(i,j).*e(i,j); % sets nu=0 outside ellipse
    
    ou=u;
    u=nu;
    t=t+dt;
    surf(x,y,u)
    colormap('jet');
    title(['t=',num2str(t,2)])
    drawnow
end
    xlabel('x')
    ylabel('y')
    zlabel('u')
    title(['2D wave equation,t=',num2str(t),', imax=',num2str(imax),', jmax=',num2str(jmax)])

% e=(x'/a).^2+(y'/b).^2<1;
% figure(2); contour(e) % to check boundary of elliptic domain
% figure(3); spy(e)     % to check elliptic domain
