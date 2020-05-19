function [P,beta] = fatigue_probability(mu_D, std_D,T)
%E_D: 
%T: Extrapolation period



T0=10;  % 10 years measurement period
E_log10_alpha = 12.76;

var_alpha=10;
var_e=0.14;
var_ee=0.1;

K=T0/T; 

E_D=mu_D;
CV_D=std_D/mu_D;

beta = (E_log10_alpha-log10(T/T0)-log10(E_D)) / sqrt(var_alpha+K*CV_D.^2+var_e+var_ee);
P = 1-normcdf(beta,1,1)