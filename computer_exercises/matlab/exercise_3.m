close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo
%%
load('../Computer Exercise Data/testdata.mat')
figure();
hold on
plot(F,X,'bx')

%%
N = length(F);
b=regress(X, [ones(N,1) F]);
xest=b(1)+b(2)*F;

plot(F,xest,'r.-');

%%
% 4.1Goodness-of-fit (ANOVA)

xest = xest;

SSR=sum((xest-mean(X)).^2);
SSE=sum((X-xest).^2);
SST=sum((X-mean(X)).^2);
% you could check if SST=SSR+SSE
if ((SST-(SSR+SSE))/(SSR+SSE)>0.01)
    warning('SST=SSR+SSE Failed');
end

MSR = SSR/1;
MSE = SSE/(N-2);

F_ratio=(MSR)/(MSE);
% P=1-fcdf(F_ratio, 1,N-2);
P=finv(0.95,1,N-2);


if F_ratio > P
    disp('X=b0 Rejected!!!');
else
    disp('X=b0 Accepted');
end
    
   
