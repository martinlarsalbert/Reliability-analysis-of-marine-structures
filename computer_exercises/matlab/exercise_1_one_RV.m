close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo

N=10000;
MU=10;
SIGMA=30;
x=normrnd(MU,SIGMA,1,N); % generate a sample of a normal R.V with N(10,302)

figure();
plot(sort(x),'.');


%%
%Compute the empirical distribution of X, F1(x), and plot it in a figure.

F1(:,1)=sort(x);
F1(:,2)=(1:N)/N;
figure();
plot(F1(:,1),F1(:,2),'b.')
hold on

%%
%Alternatively, the calculation could be done by Matlab wafo toolbox command “empdistr”,

[F X]=ecdf(x);
F2 = [X F];
plot(F2(1:50:end,1),F2(1:50:end,2),"ko");

%%
% Compute the values of mean and variance by direct formulas to get (EX1, VX1) 
% and by Matlab commands to get (EX2, VX2). Compare with the real values in the 
% simulation i.e. 

EX1=sum(x)/10000
VX1=1/(10000-1)*sum((x-EX1).^2)
EX2=mean(x)
VX2=var(x)

%%
% 1.3	Quantiles of X
%The concept of quantile is important in reliability analysis and extreme predictions. The quantile of a distribution could be defined in different ways. For the random variable X, the quantile as a number of x is defined as,
% F(x_a)=P(X<x_a)=a

x005=norminv(0.05, 10, 30)
x095=norminv(0.95, 10, 30)

%%
%2	Distributions of two or more random variables
% Let X, Y denote two independent random variables, 
% while X, Y in [0, +inf). 
% The joint pdf and cdf are represented by fX,Y(x,y) and FX,Y(x,y), respectively. 



X1=normrnd(0,1,1000,1);
X2=normrnd(0,1,1000,1);
X=X1+0.3*X2;
Y=-0.4*X1+0.1*X2;
Z=X1-X2;
R=corrcoef([X Y Z])
SIGMA=cov([X,Y,Z])





