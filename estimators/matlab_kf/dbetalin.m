function db = dbetalin(dt,u,delta,Cf,Cr,lf,lr,m,beta,r)
% This function computes dotbeta from the discretised linear bicycle model

db = -(1-dt*(Cf+Cr)/(m*u))*beta + dt*((Cf*lf-Cr*lr)/(m*u^2))*r - ...
       dt*(Cf*delta/(m*u));

end