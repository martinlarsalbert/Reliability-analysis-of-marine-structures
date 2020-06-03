close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo
%%
% Q9
%Oct.10 9 8 9.5 10 11 10.5 7 4 12.5
%Nov.13 12 14 12 9 11 6 4 2 13
%Dec. 14 13.5 13 15 14.5 16 17 20 15 16.5

wind_data = [10 9 8 9.5 10 11 10.5 7 4 12.5;
13 12 14 12 9 11 6 4 2 13;
14 13.5 13 15 14.5 16 17 20 15 16.5
];

length(wind_data)
average_wind = mean(mean(wind_data))
sum(sum(wind_data))/30

Power_wrong = 20*average_wind^3+30  % This is not correct
Power = sum(sum(20*wind_data.^3+30))/30 %...this is

%%
% Q10
S_mid = [88 130 84 78 119 73 114 109 97 93];
S_aft = [113 106 94 86 99 87 118 100 94 94];
figure(1);
plot(S_mid,S_aft,'o');
X=S_mid;
Y=S_aft;

%C(X,Y) = E(XY) - E(X)*E(Y)
mu_X = mean(X);
mu_Y = mean(Y);
C = mean(X.*Y) - mean(X)*mean(Y)
C = mean((X-mean(X)).*(Y-mean(Y))) 

C/(std(X)*std(Y))

cov(X,Y,1)
corrcoef(X,Y)

A = [50 18;
     0  24];
 Z = [1.3847 -0.6442 1.8710 0.3519 1.0569;
      0.4740 0.3073 1.2019 2.5546 0.1669]';
 
X_ = [mu_X;mu_Y] + A'*Z'
Ds = (X_).^3/(10^8);
D = sum(Ds,2)