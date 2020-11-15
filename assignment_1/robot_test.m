clearvars;clc;
Robot=load('AssignmentPredictiveMaintenance.mat');
speed=Robot.RobotArms(4).speed;
time=Robot.RobotArms(4).settling_time;
invalid_speed=speed~=4.2;
time(invalid_speed)=[];
N=length(time);
mu=1.7;
sigma=0.5;
half_alpha=0.01/2;
% sample mean
s_mean=sum(time)/length(time);
% sample variance
s_var=sum((time-s_mean).^2)/(length(time)-1);
s_sigma=sqrt(s_var);
s_z = sqrt(N).*(s_mean-mu)/s_sigma;
% p-value
pval=cdf("Normal",-abs(s_z),0,1)+(1-cdf("Normal",abs(s_z),0,1));
% convert 1-alph to x-axis 
z_half=abs(icdf('Normal',half_alpha,0,1));
% confidenc interval
ci_l=s_mean-z_half.*s_sigma/sqrt(N);
ci_r=s_mean+z_half.*s_sigma/sqrt(N);

% p = normspec([-abs(s_z),abs(s_z)],0,1,'outside');
% p = [0.025,0.975];
% ci2 = icdf('Normal',0.025,0,1);
% normspec(pval,0,1,'outside');
% [h,pvalue,ci] = ztest(v_time,mu,sigma);
% y=exp(-0.5*((x-0)/1).^2)/(1*sqrt(2*pi));

   
