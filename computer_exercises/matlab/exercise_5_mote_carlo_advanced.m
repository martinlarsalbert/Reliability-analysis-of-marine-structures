close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo
%%
% Sum of normal distributions --> normal distribution:

n=10000;
X1=normrnd(0,1,n,1);
X2=normrnd(0,1,n,1); 
Y=X1+X2;
figure();
subplot(3,1,1);
wnormplot(X1);
subplot(3,1,2);
wnormplot(X2);
subplot(3,1,3);
wnormplot(Y);

% Sum of weibul distributions is NOT weibull distribution:
lamda=1;
k=0.5;
X1=wweibrnd(lamda,k,1,n);
X2=wweibrnd(lamda,k,1,n); 
Y=X1+X2;
figure();
subplot(3,1,1);
wweibplot(X1);
subplot(3,1,2);
wweibplot(X2);
subplot(3,1,3);
wweibplot(Y);

%%
% 2.1	Correlated Normal random variables

sigma=[1.0 0.5 0.5; 0.5 2.0 0.3; 0.5 0.3 1.5];  % Covariance matrix
size=1000;
C1=chol(sigma);  %Calculate transformation matrix with Cholesky decomposition
[V D]=eig(sigma);
C2=(V*sqrt(D))';  %Calculate transformation matrix with Eigen decomposition
X1=normrnd(0,1,1,size);
X2=normrnd(0,1,1,size);
X3=normrnd(0,1,1,size);
X=[X1; X2; X3];
Z1=C1'*X;
Z2=C2'*X;

% compare the covariance of simulated R.V.s by the two methods
cov(Z1')
cov(Z2')

% Increase size to check convergence problem
size=10000;
X1=normrnd(0,1,1,size);
X2=normrnd(0,1,1,size);
X3=normrnd(0,1,1,size);
X=[X1; X2; X3];
Z1=C1'*X;
Z2=C2'*X;

% compare the covariance of simulated R.V.s by the two methods
cov(Z1')
cov(Z2')

%%
% 2.2	Non-zero mean Gaussian random vector

size=1000;
X1=normrnd(0,1,1,size);
X2=normrnd(0,1,1, size);
X3=normrnd(0,1,1,size);
X=[X1; X2; X3];
Z=C1'*X;
nominal_mean_values=[2 ; 3; 4];
Z=Z+nominal_mean_values*ones(1,size);
cov(Z')
simulated_mean_values=mean(Z,2)
nominal_mean_values

%%
%Simulating with lognormal distribution
% 1, Zero means of Z
size=1000;
X1=normrnd(0,1,1,size);
X2=normrnd(0,1,1, size);
X3=normrnd(0,1,1,size);
X=[X1; X2; X3];
Z=C1'*X;
Y=exp(Z);
cov(Y')
mean_0_exp = mean(Y',1)

%2, Z has the mean values of [2; 3; 4]
Z=Z+nominal_mean_values*ones(1,size);
Y=exp(Z);
cov(Y')
nominal_mean_values_exp = mean(Y',1)

%%
% 4.1	Expected and variance of random variables functions
%Now suppose that Z1, Z2 are correlated normal random variables. 
% The covariance matrix is sigma=[3.5 1.3; 1.3; 5], 
%and the mean values are E[Z1, Z2]=[1.5  5]. 
%Compute the mean and variance of the function f(Z1, Z2)=Z1^2+Z1*Z2+Z2^2.

sigma=[3.5 1.3; 1.3 5];
mea=[1.5; 5];
size=10000;
C=chol(sigma);
X1=normrnd(0,1,1,size);
X2=normrnd(0,1,1,size);
X=[X1; X2];
Z=C'*X;
Z=Z+mea*ones(1,size);
Z1=Z(1,:); 
Z2=Z(2,:);
f=Z1.^2+Z1.*Z2+Z2.^2;
[mean(f) var(f)]

%%
% 4.2	Confidence interval of (function of) correlated random variables
% If a random variable is normally distributed, the confidence interval can be computed using only the mean and variance of the random variable. 
%However, when the random variable is not normally distributed, or even as above a complex function of several correlated normally distributed random variables, 
%the computation of confidence interval should be computed using the original definition, 
%i.e. P(x_α/2<X<x_1-α/2)=1-α. For example, the 95% confidence interval of
%the above function of two correlated random variables is:
f=sort(f);
x025=f(size*0.025);
x975=f(size*0.975);
CI1=[x025,  x975]

figure()
wnormplot(f);
CI2=[mean(f)-1.96*std(f)  mean(f)+1.96*std(f)]  % f is clearly not normal distributed and calculating the convidence intervalls this way will not be correct.








