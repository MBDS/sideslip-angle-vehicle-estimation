%% Beta estimation from Stanford data
% Author: A. Leanza, 04/2021

close all
clear
clc

addpath C:\Users\Antonio\Documents\CartelleCorrenti\Functions

rng default



%% Load data

% Experimental data
% Car data
lf = 1.33;  % m
lr = 1.07;  % m
m = 982;  % kg
a = lf+0.841;  b = lr+0.849;  % m
h = 1.115;  % m
track = 1.35;  % m  (che poi sarebbe tf e/o tr)
Len = 4.09;  % m
% Jz = m/3*(a^2-a*b+b^2)  % kgm^2  ANTONIO
Jz = m/12*(1.7^2+4.09^2);  % kgm^2  LETTERATURA
% Jz = m*lf*lr  % kgm^2  LETTERATURA
% Jz = 0.17*h^(1/3)*Len^1.46  % kgm^2  LETTERATURA  
% Jz = m*track*2.4/2.2048;  % kgm^2  LETTERATURA  
% Jz = 1500;
L = lr + lf;  % m
tf = 1.35;  % m
tr = tf;  % m
% Measures        
race = 1;
switch race
    case 1
        load('20140222_01_01_03_250lm.mat')
        start = 15000;  fin = 70000;
    case 2
        load('20140222_02_01_03_250lm.mat')
        start = 25000;  fin = 75000;
end
dt = 0.01;  % 100 Hz
t = insData.ayCG.time;  % s
ay = insData.ayCG.value;  % m/s^2
ayFilt = insData.ayCGFilt.value;  % m/s^2 already filtered (interesting!)
ax = insData.axCG.value;  % m/s^2
axFilt = insData.axCGFilt.value;  % m/s^2 already filtered (interesting!)
r = insData.yawRate.value*pi/180;  % rad/s
rFilt = insData.yawRateFilt.value*pi/180;  % rad/s already filtered (interesting!)
dr = insData.yawAngAcc.value*pi/180;  % rad/s^2
drFilt = insData.yawAngAccFilt.value*pi/180;  % rad/s^2 already filtered (interesting!)
tau = 13.529/180;
% delta = driverData.handwheelAngle.value*tau*pi/180;  % rad with 1kHz
delta = tireData.roadWheelAngle.value*pi/180;  % rad with 1kHz
delta = delta(1:10:end);  % rad from 1kHz to 100 Hz

delta_r = tireData.roadWheelAngle.value*pi/180;  % delta road
delta_r = delta_r(1:10:end);

deltaFL = tireData.roadWheelAngleFL.value*pi/180;
deltaFL = deltaFL(1:10:end);
deltaFR = tireData.roadWheelAngleFR.value*pi/180;
deltaFR = deltaFR(1:10:end);

% Rfl = 0.2976;  % m
% Rfr = 0.2976;  % m
% Rrl = 0.3199;  % m
% Rrr = 0.3199;  % m

vx = insData.vxCG.value;
vy = insData.vyCG.value;

beta_true = insData.sideSlip.value*pi/180;

wsFL = tireData.wheelSpeedFL.value;
wsFR = tireData.wheelSpeedFR.value;
wsRL = tireData.wheelSpeedRL.value;
wsRR = tireData.wheelSpeedRR.value;

% Elimino il tratto a vx = 0.
% start = 4452;  fin = 32216;
t = t(start:fin);
ay = ay(start:fin);
ayFilt = ayFilt(start:fin);
ax = ax(start:fin);
axFilt = axFilt(start:fin);
r = r(start:fin);
rFilt = rFilt(start:fin);
delta = delta(start:fin);
vx = vx(start:fin);
vy = vy(start:fin);
beta_true = beta_true(start:fin);
wsFL = wsFL(start:fin);
wsFR = wsFR(start:fin);
wsRL = wsRL(start:fin);
wsRR = wsRR(start:fin);


% plot(t,wsFL,t,wsFR,t,wsRL,t,wsRR,t,vx)

if 0
% Sensors noise
% Design the high pass filter
d = fdesign.highpass('Fst,Fp,Ast,Ap',0.15,0.25,60,1);
% designmethods(d)
Hd = design(d,'equiripple');
% fvtool(Hd)
Noise_wsFL = filter(Hd,wsFL);
Noise_wsFR = filter(Hd,wsFR);
Noise_wsRL = filter(Hd,wsRL);
Noise_wsRR = filter(Hd,wsRR);


subplot(4,1,1)
plot(t,Noise_wsFL)
subplot(4,1,2)
plot(t,Noise_wsFR)
subplot(4,1,3)
plot(t,Noise_wsRL)
subplot(4,1,4)
plot(t,Noise_wsRR)
end


T = length(t);  


%% KF for lateral beta estimation

% Sensors noise
% Design the high pass filter
d = fdesign.highpass('Fst,Fp,Ast,Ap',0.15,0.25,60,1);
% designmethods(d)
Hd = design(d,'equiripple');
% fvtool(Hd)
Noise_ay = filter(Hd,ay);
s_ay = std(Noise_ay);
% s_ay = 2.5; 
Noise_r = filter(Hd,r);
s_r = std(Noise_r);  % rad/s
Noise_delta = filter(Hd,delta);
s_d = std(Noise_delta)*1000;  % rad

 
Cf = 0.7*10^5;
Cr = 1.2*10^5;


% Linear KF
Q = s_d^2;
R = diag([s_ay^2 s_r^2]);
X = [0; 0];
P = diag([10^4 10^4]);
Y = [ay'; r'];

for i =2:T
    
    Ad = [1-dt*(Cf+Cr)/(m*vx(i-1)) -dt*(1+(lf*Cf-lr*Cr)/(m*vx(i-1)^2));
          -dt*(lf*Cf-lr*Cr)/Jz 1-dt*(lf^2*Cf+lr^2*Cr)/(Jz*vx(i-1))]; 
    
    Bd = dt*[Cf/(m*vx(i-1)); lf*Cf/Jz]; 
    
    Gd = Bd/dt; 
    
    Qd = Bd*Q*Bd'; 
    
%     Cd = [-(Cf+Cr)/(m*vx(i)) -(lf*Cf-lr*Cr)/(m*vx(i)^2);
%           0 1];
Cd = [-(Cf+Cr)/m -(lf*Cf-lr*Cr)/(m*vx(i));
          0 1];
      
%     Dd = [Cf/(m*vx(i)); 0];  
Dd = [Cf/m; 0];
    
    Rd = R;
    
    [X,P] = LinearKF(Ad,Bd,Cd,Dd,delta(i-1),delta(i),Y(:,i),X,P,Qd,Rd);
    
    States(:,i) = X;
      
end

beta = States(1,:); beta = beta';

RMSEbeta = sqrt(mean((beta-beta_true).^2))*180/pi;
RMSEr = sqrt(mean((States(2,:)'-r).^2))*180/pi;

figure(2)
set(2,'units','normalized','position',[0.0012 0.2142 0.7713 0.7150],....
     'name','Fig1','numbertitle','off')
subplot(2,1,1) 
plot(t,beta_true*180/pi,'k','linewidth',2), hold on
plot(t,beta*180/pi,'color',[0.5 0.5 0.5],'linewidth',2)
text(t(1)*1.1,-3.5,strcat('RMSE=',num2str(RMSEbeta,'%.2f\n'),' deg'),'Fontsize',16);
xlabel('t [s]'), ylabel('\beta [deg]')
legend('Actual','Estimated')
axis([t(1) t(end) -5 5])
set(gca,'FontSize',14)
subplot(2,1,2)
plot(t,r*180/pi,'k','linewidth',2), hold on
plot(t,States(2,:)*180/pi,'color',[0.5 0.5 0.5],'linewidth',2)
text(t(1)*1.1,30,strcat('RMSE=',num2str(RMSEr,'%.2f\n'),' deg/s'),'Fontsize',16);
xlabel('t [s]'), ylabel('r [deg/s]')
legend('Actual','Estimated')
xlim([t(1) t(end)])
set(gca,'FontSize',14)


