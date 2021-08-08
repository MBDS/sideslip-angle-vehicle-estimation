function dr = dyawlin(dt,u,delta,Cf,Cr,lf,lr,Jz,beta,r)
% This function computes dotyaw from the discretised linear bicycle model

dr = -(1-dt*(Cf*lf^2+Cr*lr^2)/(Jz*u))*r + ...
     dt*(Cf+lf-Cr*lr)/Jz*beta - dt*(Cf*lf*delta)/Jz;

end