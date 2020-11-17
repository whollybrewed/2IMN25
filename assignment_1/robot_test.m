clearvars;clc;
Robot=load('AssignmentPredictiveMaintenance.mat');
for index=1:28
    speed=Robot.RobotArms(index).speed;
    time=Robot.RobotArms(index).settling_time;
    invalid_speed=speed~=4.2;
    time(invalid_speed)=[];
    N=length(time);
    mu=1.7;
    sigma=0.5;
    sig_lv=0.01;
    % sample mean
    s_mean=sum(time)/length(time);
    results(index).s_mean=s_mean;
    % sample variance
    s_var=sum((time-s_mean).^2)/(length(time)-1);
    s_sigma=sqrt(s_var);
    results(index).s_sigma=s_sigma;
    s_z = sqrt(N).*(s_mean-mu)/s_sigma;
    % p-value
    pval=cdf("Normal",-abs(s_z),0,1)+(1-cdf("Normal",abs(s_z),0,1));
    results(index).pval=pval;
    if pval < sig_lv
        report(index).ztest=1;
    else
        report(index).ztest=0;
    end
    % half width  
    half_width=abs(icdf('Normal',sig_lv/2,0,1)).*s_sigma/sqrt(N);
    % confidenc interval
    ci_l=s_mean-half_width;
    ci_r=s_mean+half_width;
    results(index).ci=[ci_l ci_r];
    %chi squared testing
    k = 100;
    x = norminv(1/k:1/k:(k-1)/k,mu,sigma);
    pd = makedist('Normal','mu',mu,'sigma',sigma);
    y=diff(cdf(pd,x));
    y=[(1-sum(y))/2 y (1-sum(y))/2];
    expcount=N*y;
    edges=[-inf x inf];
    bincount=histcounts(time,edges);
    chisquare=sum(((bincount-expcount).^2)./expcount);
    c_val_chi=icdf("Chisquare",0.99,k-1)
    results(index).chisquare=chisquare;
    if chisquare > c_val_chi
        report(index).chi2=1;
    else
        report(index).chi2=0;
    end
    %chi-squared test for variance
    var_chisquare=(N-1).*s_var/(sigma.^2);
    c_val_chi_var=icdf("Chisquare",0.99,N-1)
    results(index).var_chisquare=var_chisquare;
    if var_chisquare > c_val_chi_var
        report(index).chi2_var=1;
    else
        report(index).chi2_var=0;
    end
    % KS test
    quantile=sort(time.');
    emp_cdf=[];
    for n=1:N
        emp_cdf=[emp_cdf sum(quantile(:)<=quantile(n))/N];
    end
    norm_cdf=cdf(pd,quantile);
    d_plus=max(emp_cdf-norm_cdf);
    d_minus=max(norm_cdf-(emp_cdf-1/N));
    d_n=max(d_plus,d_minus);
    c_val_ks=1.628/(sqrt(N)+0.12+(0.11/sqrt(N)))
    results(index).d_n=d_n;
    if d_n > c_val_ks
        report(index).ks=1;
    else
        report(index).ks=0;
    end
end

   
