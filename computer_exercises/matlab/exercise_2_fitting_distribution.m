close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo

%%
% 1	Plot the data to check which type of distributions
% Assume that we have a set of observations x1, x2, … , xn. Before we estimate any parameters, we must convince ourselves that the observations originate from the right family of distributions, e.g. normal, Gumbel, or Weibull. 
% One way to get a rough idea, about which family of distributions may be suitable, is to display the observations in a probability plot. If you suspect that the data originate from, for instance, a normal distribution, then you should make a normal probability plot; if you instead suspect a Gumbel distribution, and then make a Gumbel probability plot. If, in the plot, the observations seem to line up well along a straight line, it indicates that the chosen distribution for the probability plot indeed might serve as a good model for the observations. 
% Statistics Toolbox provides normplot (for normal distribution), weibplot (for Weibull distribution); the WAFO toolbox furnishes you with wgumbplot (for Gumbel distribution). Acquaint yourself with the above-mentioned commands, for example

dat1=randn(2000,1); % generate normal random variable, referred later on.
figure();
normplot(dat1)

figure();
wweibplot(dat1) % Any error-message?

%%
dat2=rand(3000,1); % generate Uniform distributed random variable.
figure();
normplot(dat2)
figure();
wgumbplot(dat2)

dat3=wweibrnd(2,2.3,1,3000); % generate Weibull distributed random variable.
figure();
wweibplot(dat3)
figure();
wgumbplot(dat3)

%%
% 1.1	Measurements of significant wave height Hs in the Atlantic Ocean
% In oceanography and marine technology, statistical extreme-value theory has been used to a great extent. In design of offshore structures knowledge about “extreme” condition is important.
% In the numerical examples above, we used artificial data, simulated from a distribution which we could control. We will now consider real measurements from the Atlantic Ocean. The data set contains so-called significant wave heights. 
atl=load('../Computer Exercise Data/atlantic.dat');
figure();
plot(atl,'.')

figure()
subplot(4,1,1);
normplot(atl)

subplot(4,1,2);
normplot(log(atl))

subplot(4,1,3);
wgumbplot(atl)

subplot(4,1,4);
wweibplot(atl)

%%
% 2.2	Estimation errors (Mean and variance of the estimated parameters)
% In the WAFO toolbox the ML-method has been implemented in “wgumbfit”, “wweibfit”, and “wraylfit” for the purpose of estimating the parameters in a Gumbel, Weibull, and Rayleigh distributions, respectively.
% •	Please compute the two parameters by the WAFO routine “wgumbfit”, and then compare with that got from above direct calculation.

figure()
[phat, covm]=wgumbfit(atl);
