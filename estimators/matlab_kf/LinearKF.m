function [X,P] = LinearKF(Ad,Bd,Cd,Dd,u_,u,Y,X_,P_,Qd,Rd)
% This function computes the Linear Kalman Filter

% Prediction
Xpred = Ad*X_ + Bd*u_;
Ppred = Ad*P_*Ad' + Qd;

% Correction
S = Cd*Ppred*Cd' + Rd;
K = Ppred*Cd'*S^(-1);
X = Xpred + K*(Y-Cd*Xpred-Dd*u);
P = (eye(length(X))-K*Cd)*Ppred;