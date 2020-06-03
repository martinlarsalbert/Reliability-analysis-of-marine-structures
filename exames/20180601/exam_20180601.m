close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo
%%
% Q8
% Hasofer & Lind:

a_start = -0.58*[1,1,1];
beta_start = 3;

a = a_start;
beta = beta_start;

as = [];
betas = [];
ks = [];

for i=1:10
    [a, beta, k] = hasofer(a,beta);
    as(i,:) = a;
    betas(i)=beta;
    ks(i)=k;
end;
a
beta
k

beta_HF = beta

figure(1);
plot(ks);
legend('k');

figure(2);
plot(beta);
legend('beta');

figure(3);
plot(as);
legend('a_1','a_2','a_3');
 
Pf = 1 - normcdf(beta_HF,0)
