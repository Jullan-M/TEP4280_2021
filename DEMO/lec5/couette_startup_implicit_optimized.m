% Implicit Euler method with central FDM (BTCS) for startup of Couette flow

clear all
close all
clc

jmax=21;
y=linspace(0,1,jmax);
dy=1/(jmax-1);

r=1e100; %von Neumann number
dt=r*dy^2;
%tstop=1; % end time
tstop= dt; % fzero(@u_dev,[0 1])
nmax=ceil(tstop/dt);
dt=tstop/nmax;
r=dt/dy^2;

% Initial condition
u=zeros(jmax,1);
u(jmax)=1;
unew=u;

% A=(1+2*r)*diag(ones(jmax,1),0)-r*diag(ones(jmax-1,1),-1)-r*diag(ones(jmax-1,1),1);
% A(1,1)=1; A(1,2)=0;
% A(jmax,jmax-1)=0; A(jmax,jmax)=1;

A=(1+2*r)*diag(ones(jmax-2,1),0)-r*diag(ones(jmax-3,1),-1)-r*diag(ones(jmax-3,1),1);
    
% Time loop for FTCS
for n=1:nmax
    % BTCS enhanced
    b = u(2:jmax-1);
    b(end) = b(end)+r;
    u(2:jmax-1) =A\b;
    t=n*dt;

    u2=zeros(jmax,1);
    for l=1:15
       u2=u2-sin(l*pi*(1-y'))*exp(-(l*pi)^2*t)/l;
    end
    u2=y'+u2*2/pi;

    u3=zeros(jmax,1);
    fac=0.5/sqrt(t);
    for l=0:5
       u3=u3+erfc((2*l+1-y')*fac)-erfc((2*l+1+y')*fac);
    end

    plot(u,y,'x',u2,y,'g',u3,y,'r')
    xlim([0 1])
    title(['t = ',num2str(t,2)])
    drawnow
end
xlabel('u(y)')
ylabel('y')
hold on
u_ex=y';
plot(u_ex,y,'--c')
legend('BTCS','exact (2)','exact (3)','exact steady state solution','Location','best')
hold off
