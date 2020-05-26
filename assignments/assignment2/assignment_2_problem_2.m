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


N=30;
result_course_1 = simulate(C1, data_course_1.MEAN, N);
result_course_2 = simulate(C2, data_course_2.MEAN, N);

%%
figure(3);
p1 = plot(result_course_1.Hs);
for i=1:length(p1)
    p1(i).Color(4) = 0.25;
end
hold on;

plot(exp(data_course_1.MEAN),'r-','LineWidth',2);
legend('Mean significant wave height'); 

xlabel('Sea state number');
ylabel('Significant wave height [m]');
title('Significant wave height Monte Carlo Simulations Course 1 ');


figure(4);
[X,Y] = meshgrid(1:length(data_course_1.MEAN),1:length(data_course_1.MEAN));
contourf(X,Y,data_course_1.SIGMA);
title('Covariance matrix Course 1');
xlabel('Sea state number');
ylabel('Sea state number');
colorbar();

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

%%
figure(5);
normplot(log(D_voy_));
y_=log(D_voy_);

figure(6);
histogram(D_voy_,100,'Normalization','probability');

figure(1);
plot(ns,means,'o');
ylims = ylim();
ylim([0,ylims(2)]);
grid;
xlabel('N');
ylabel('mean(D_{voy})');
title('Check mean value convergence');

figure(2);
plot(ns,vars,'o');
ylims = ylim();
ylim([0,ylims(2)]);
grid;
xlabel('N');
ylabel('var(D_{voy})');
title('Check variance convergence');

%% ToDo: Check confidence, how?

%%
% Fatigue reliability analysis
E_D = means(end);

mu_D = means(end);
var_D = vars(end);
std_D = sqrt(var_D);

[P_5,beta_5] = fatigue_probability(mu_D, std_D, 5)
[P_10,beta_10] = fatigue_probability(mu_D, std_D, 10)
[P_100,beta_100] = fatigue_probability(mu_D, std_D, 100)
[P_1000,beta_1000] = fatigue_probability(mu_D, std_D, 1000)


%%
%Save figures:
exportgraphics(figure(3),'mean_wave_course_1.pdf');
exportgraphics(figure(4),'covariance_course_1.pdf');







