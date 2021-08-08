function [] = compare_results(filesPrefix)
% Usage:  compare_results('fg_Win5_');
%
dir='../datasets/20140222_01_01_03_250lm';

dt = 0.01;  % 100 Hz
t = load(sprintf('%s/t.txt',dir));
beta_true = load(sprintf('%s/beta_true.txt',dir));
r = load(sprintf('%s/yawRate.txt',dir));
vx = load(sprintf('%s/vx.txt',dir));


%% Load c++ results
t_beta=load(sprintf('%sestimated_beta.txt',filesPrefix));
beta=t_beta(:,2);

t_yawrate=load(sprintf('%sestimated_yawrate.txt',filesPrefix));
yawrate=t_yawrate(:,2);

N = size(t_beta,1);

%t_yawrate_obs=load('build-Release/out_obs_yawrate.txt');
%yawrate_obs=t_yawrate_obs(:,2);


%% Plots 

RMSEbeta = sqrt(mean((beta(1:N)-beta_true(1:N)).^2))*180/pi;
RMSEr = sqrt(mean((yawrate(1:N)-r(1:N)).^2))*180/pi;

afigure()
%set(2,'units','normalized','position',[0.0012 0.2142 0.7713 0.7150],....     'name','Fig1','numbertitle','off')
subplot(2,1,1) 
plot(t(1:N),beta_true(1:N)*180/pi,'k','linewidth',2), hold on
plot(t(1:N),beta(1:N)*180/pi,'color',[0.5 0.5 0.5],'linewidth',2)
text(t(1)*1.1,-3.5,strcat('RMSE=',num2str(RMSEbeta,'%.2f\n'),' deg'));
xlabel('t [s]'), ylabel('\beta [deg]')
legend('Actual','Estimated')
axis([t(1) t(N) -5 5])

subplot(2,1,2)
plot(t(1:N),r(1:N)*180/pi,'k','linewidth',2), hold on
plot(t(1:N),yawrate(1:N)*180/pi,'color',[0.5 0.5 0.5],'linewidth',2)
text(t(1)*1.1,30,strcat('RMSE=',num2str(RMSEr,'%.2f\n'),' deg/s'),'Fontsize',16);
xlabel('t [s]'), ylabel('r [deg/s]')
legend('Actual','Estimated')
xlim([t(1) t(N)])

%subplot(3,1,3)
%plot(t(1:N),vx(1:N),'k','linewidth',2);
%title('vx');
%axis tight;

end
