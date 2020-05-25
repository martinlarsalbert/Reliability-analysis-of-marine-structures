close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo

%%
% 1  Fitting distribution to yearly maxima of Hs
% The measurement  data is recorded from   1979-01-01  to  2011-12-31,
% measuring once per day at 0:00:00.
data = load('Hs_79_2011.mat');
Hs = data.Hs;

figure(1);
plot(Hs);

% Compute and plot the empirical distribution of the yearly maxima:
number_of_years = 2012-1979;
year_length = (length(Hs) - mod(length(Hs), number_of_years))/number_of_years;

year_start_index=1:year_length:length(Hs);
year_stop_index=year_start_index(2:end)-1;
year_start_index=year_start_index(1:end-1);

years = {};
yearly_maximas = [];
for i=1:length(year_start_index)
    year = struct();
    start=year_start_index(i);
    stop=year_stop_index(i);
    year.Hs = Hs(start:stop);
    yearly_maximas(end+1)=max(year.Hs);
    years{i}=year;
end
yearly_maximas=sort(yearly_maximas);

figure(2);
for i=1:length(years)
    year=years{i};   
    [x,y] = empirical(year.Hs);
    p1=plot(x,y,'b-');
    hold on;
    p1.Color(4) = 0.20;
end
title('Empirical yearly distibution of H_s 1979-2012');
xlabel('x=H_s [m]');
ylabel('P(H_s < x)');

figure(3);
[x,y] = empirical(yearly_maximas);
plot(x,y,'b-');
title('Empirical distribution of the yearly maxima');
xlabel('x_T=max(H_s)');
ylabel('P(H_s < x_t)');

%%
%ToDo: Calculate ML Estimates.

%%
figure(4);
[phat, cov] = wgumbfit(yearly_maximas);  
log_likelihood_gumb = sum(log(wgumbpdf(yearly_maximas,phat(1),phat(2))));

figure(5);
method='ML';
[phat,cov] = wgevfit(yearly_maximas,method);
log_likelihood_gev = sum(log(wgevpdf(yearly_maximas,phat(1),phat(2),phat(3))));
DEV = 2*(log_likelihood_gev - log_likelihood_gumb);
Xi2 = wchi2inv(0.95,1);
if DEV > Xi2
    disp('Gumbel is rejected');
else
    disp('Gumbel is accepted');
end;

%%
% 2  Extreme prediction of 1000-year wave
[phat, cov] = wgumbfit(yearly_maximas,0);  
exceedance_probability = 1-wgumbcdf(yearly_maximas,phat(1),phat(2));
figure(6);
[x,y] = empirical(yearly_maximas);
plot(x,1-y,'x');
hold on;

plot(yearly_maximas,exceedance_probability);
ylabel('Exceedance probability');
xlabel('H_s [m]');
legend('Empirical distribution','Fitted gumbel distribution');

lambda_t = 1000;
P=1-1/lambda_t;
H_s_1000=wgumbinv(P,phat(1),phat(2))

%%
% 3  Uncertainty analysis
sigma = phat(1);
mu = phat(2);
H_s_1000_ = mu + sigma*log(lambda_t)
n=length(years);
sigma_epsilon = sqrt(1.11*sigma^2/n + log(lambda_t)^2*0.61*sigma^2/n + 0.52*log(lambda_t)*sigma^2/n) 
alpha = 0.05;
H_s_intervall = [H_s_1000_-normcdf(alpha/2,0,1)*sigma_epsilon, H_s_1000_+normcdf(1-alpha/2,0,1)*sigma_epsilon]

%%
%Save figures:
exportgraphics(figure(2),'empirical_years.pdf');
exportgraphics(figure(3),'empirical_yearly_maxima.pdf');
exportgraphics(figure(4),'gumbel_fit.pdf');
exportgraphics(figure(5),'GEV_fit.pdf');
exportgraphics(figure(6),'exceedance_probability.pdf');

