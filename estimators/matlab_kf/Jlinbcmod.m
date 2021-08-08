function Ma = Jlinbcmod(dt,u,Cf,Cr,lf,lr,m,Jz)
% This function computes the jacobian of the discretised linear bicycle model

Ma = [-(1-dt*(Cf+Cr)/(m*u)) dt*((Cf*lf-Cr*lr)/(m*u^2)+1) 1 0;
      dt*(Cf*lf-Cr*lr)/(Jz) -(1-dt*(Cf*lf^2+Cr*lr^2)/(Jz*u)) 0 1;
      0 0 0 -1 ;
      0 0 (Cf+Cr)/m (Cf*lf-Cr*lr)/(m*u)];  
  
end
