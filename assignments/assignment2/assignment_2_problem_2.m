close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo

%%
data_course_1 = load('COURSE1.mat');  % mean values and covarianve matrix for course 1
[V D]=eig(data_course_1.SIGMA);
C1=(V*sqrt(D))';  %Calculate transformation matrix with Eigen decomposition

data_course_2 = load('COURSE2.mat');  % mean values and covarianve matrix for course 1
[V D]=eig(data_course_2.SIGMA);
C2=(V*sqrt(D))';  %Calculate transformation matrix with Eigen decomposition


N=100;
result_course_1 = simulate(C1, data_course_1.MEAN, N);
result_course_2 = simulate(C2, data_course_2.MEAN, N);
%%
D_voy = result_course_1.D_voy + result_course_2.D_voy; 
mean(D_voy)
var(D_voy)

%%
% Convergence...

means = [];
vars = [];

ns = [100,1000,3000,5000,7000,10000,20000];
for i=1:length(ns)
    n=ns(i);
    r1 = simulate(C1, data_course_1.MEAN, n);
    r2 = simulate(C2, data_course_2.MEAN, n);
    D_voy_ = r1.D_voy + r2.D_voy; 
    means(i)=mean(D_voy_);
    vars(i)=var(D_voy_);
     
end

figure(1);
plot(ns,means,'o');
ylims = ylim();
ylim([0,ylims(2)]);
grid;
xlabel('N');
ylabel('mean(D_voy)');
title('Check mean value convergence');

figure(2);
plot(ns,vars,'o');
ylims = ylim();
ylim([0,ylims(2)]);
grid;
xlabel('N');
ylabel('var(D_voy)');
title('Check variance convergence');

%% ToDo: Check confidence, how?

%%
% Fatigue reliability analysis
T=5;   % Extrapolation period
T0=10;  % 10 years measurement period

E_log10_alpha = 12.76;
E_D = means(end);
var_alpha=10;
var_e=0.14;
var_ee=0.1;

K=T0/T; 
std_D=sqrt(vars(end));
mu_D=means(end);
CV_D=std_D/mu_D;

beta = (E_log10_alpha-log10(T/T0)-log10(E_D)) / sqrt(var_alpha+K*CV_D.^2+var_e+var_ee);


% Save figures:
%exportgraphics(figure(1),'data.pdf');







