close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo
%%
dat1=randn(2000,1); % generate normal random variable, referred later on.
figure(1);
normplot(dat1);

weibplot(dat1) % Any error-message?
dat2=rand(3000,1); % generate Uniform distributed random variable.
normplot(dat2)
wgumbplot(dat2)
dat3=wweibrnd(2,2.3,1,3000); % generate Weibull distributed random variable.
wweibplot(dat3)
wgumbplot(dat3)
