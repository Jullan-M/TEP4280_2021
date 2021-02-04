% FTCS for startup of Couette flow

clear all
close all
clc

jmax=21;
y=linspace(0,1,jmax);
dy=1/(jmax-1);

r=0.5; %von Neumann number
dt=r*dy^2;
%tstop=1; % end time
tstop=fzero(@u_dev,[0 1])
nmax=ceil(tstop/dt);
dt=tstop/nmax;
r=dt/dy^2;

% Initial condition
u=zeros(jmax,1);
u(jmax)=1;
unew=u;

j=[2:jmax-1];
% Time loop for FTCS
for n=1:nmax
   % FTCS
%    for j=2:jmax-1
%        unew(j)=(1-2*r)*u(j)+r*(u(j-1)+u(j+1));
%    end
%    u=unew;
   u(j)=(1-2*r)*u(j)+r*(u(j-1)+u(j+1)); % Array operations are much faster than for loops!
  
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
legend('FTCS','exact (2)','exact (3)','exact steady state solution')
hold off

   


