% Function to find time t when the exact solution of Couette startup flow 
% is equal to 1/3 at y=1/2

function d = u_dev(t)

   u3=0;
   y=0.5;
   fac=0.5/sqrt(t);
   for l=0:5
       u3=u3+erfc((2*l+1-y)*fac)-erfc((2*l+1+y)*fac);
   end
   d=u3-1/3;