close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo
%%

%%
% Question 9.1
xs = [55,64,89,120,60,75,85,62,95,90];
lambda = mean(xs);  % According to likelihood estimition of exponential distribution

% Question 9.2
m = 3;
alpha = 10^12.76;
n_cycles = 10000; % Total number of cycles during one year.

n=10;  % Divide the into 10 stress ranges
stress_min = 0;
stress_max = 200; % [MPa]
x_limits = linspace(stress_min,stress_max,n+1);
F = 1-exp(-x_limits/lambda);  % Probability of stress lower than x_limits

figure();
plot(x_limits,F,'o-');

P=diff(F);  % Probability of stress in the stress ranges.
x = x_limits(1:end-1)+diff(x_limits)/2;  %Mid stress in the ranges

Ncrit = exp(log(alpha) - m*log(x))  % Critical number of cycler per stress range
N = P*n_cycles;  % Expected number of cycles per stress range

Dn = N./Ncrit;
D = sum(Dn);


