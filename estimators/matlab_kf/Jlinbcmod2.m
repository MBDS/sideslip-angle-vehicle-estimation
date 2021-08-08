function Ma = Jlinbcmod2(dt,u1,u2,Cf,Cr,lf,lr,m,Jz)
% This function computes the jacobian of the discretised linear bicycle model
% Output rows & columns order:
%
%         [  d_eb/d_b{k-1}  d_eb/d_r{k-1}  d_eb/d_b{k}  d_eb/d_r{k} ]
%         [  d_er/d_b{k-1}  d_er/d_r{k-1}  d_er/d_b{k}  d_er/d_r{k} ]
%  Ma =   [  d_ep/d_b{k-1}  d_ep/d_r{k-1}  d_ep/d_b{k}  d_ep/d_r{k} ]
%         [  d_eay/d_b{k-1} d_eay/d_r{k-1} d_eay/d_b{k} d_eay/d_r{k} ]
%


Ma = [-(1-dt*(Cf+Cr)/(m*u1)) dt*((Cf*lf-Cr*lr)/(m*u1^2)+1) 1 0;
      dt*(Cf*lf-Cr*lr)/(Jz) -(1-dt*(Cf*lf^2+Cr*lr^2)/(Jz*u1)) 0 1;
      0 0 0 -1 ;
      0 0 (Cf+Cr)/m (Cf*lf-Cr*lr)/(m*u2)];  
  
end
