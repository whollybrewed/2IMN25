clearvars;clc;
Robot=load('AssignmentPredictiveMaintenance.mat');
speed=Robot.RobotArms(19).speed;
time=Robot.RobotArms(19).settling_time;
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
% half width  
half_width=abs(icdf('Normal',half_alpha,0,1)).*s_sigma/sqrt(N);
% confidenc interval
ci_l=s_mean-half_width;
ci_r=s_mean+half_width;
%chi squared testing
bin = 100;
x = norminv(1/bin:1/bin:(bin-1)/bin,mu,sigma);
pd = makedist('Normal','mu',mu,'sigma',sigma);
y=diff(cdf(pd,x));
y=[(1-sum(y))/2 y (1-sum(y))/2];
expcount=N*y;
edges=[-inf x inf];
bincount=histcounts(time,edges);
chisquare=sum(((bincount-expcount).^2)./expcount);
icdf("Chisquare",0.99,bin-1);


   
