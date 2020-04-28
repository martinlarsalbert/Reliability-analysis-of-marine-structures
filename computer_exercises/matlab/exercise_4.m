close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo
%%
% First, we will use the above theory to generate 100 Weibull-distributed numbers:
x_u = rand(100,1);

lambda=2; theta=0; k=3.6;
x_weibull=theta+lambda*(-log(1-x_u)).^(1/k);
figure();
plot(x_weibull,'b.'); hold on

% Or we could generate 100 normally-distributed numbers viz.
mu=10; sigma=3;

x_normal=mu+sigma*norminv(x_u);
plot(x_normal,'bo')

% Eventually, generate 100 Gumbel-distributed numbers:
mu=3.6; beta=2;
x_gumbel=mu-beta*log(-log(x_u));
plot(x_gumbel,'r*'); hold on

legend('Weibull','Normal','Gumbel');


%%
figure();
wweibplot(x_weibull);

figure();
wnormplot(x_normal);

figure();
wgumbplot(x_gumbel);

%%
% 3.1	Approximating the value of π

%first to draw the circle and square
theta=0:0.1:360;
x=sind(theta); 
y=cosd(theta);
figure();
hold on;
plot(x,y); 
plot([-1 -1],[-1 1],'k-');
plot([-1 1],[-1 -1],'k-');
plot([1 1],[-1 1],'k-');
plot([-1 1],[1 1],'k-');
size=100;
x1=-1+2*rand(size,1);
y1=-1+2*rand(size,1);
plot(x1,y1,'b.');
inds=find(x1.^2+y1.^2<=1);
R=length(inds)/size;
PI=4*R

y=sqrt(1-x1.^2);
V=(1-(-1))^2;
PI2=mean(y)*V

%%
%first to draw the circle and square
theta=0:0.1:90;
x=sind(theta); 
y=cosd(theta);
figure();
hold on;
plot(x,y); 
size=400;
x1=rand(size,1);
y1=rand(size,1);
plot(x1,y1,'b.');
inds=find((x1.^2+y1.^2)<1);
R=length(inds)/size;
PI=4*R

y=sqrt(1-x1.^2);
V=1;
PI2=4*mean(y)*V




