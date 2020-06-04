close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo
%%
Vs= [5
6
7
8
9
10];

Ps = [650
1000
1300
2000
2500
3200];

xs=Vs.^3;
ys=Ps;
n=length(xs);

b1=(sum(xs.*ys)-sum(xs)*sum(ys)/n)/sum((xs-mean(xs)).^2);
b0=mean(ys)-b1*mean(xs);

x=9.5^3;
P = b0+b1*x

V_=linspace(min(Vs),max(Vs),20);
x_=V_.^3;
P_=b0+b1*x_;


figure(1);
plot(Vs,Ps,'o');
hold on;
plot(V_,P_);
xlabel('V');
ylabel('P');
legend('Mesurements','Empirical formula');
