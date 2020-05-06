close all;
clear all;
clc;
addpath('../wafo_old')
initwafo

%%
load('Virkler.mat');  % Loading into Virk matrix

% Note that only the first 136 rows could be used for the analysis in thisp roject. 
data = Virk(1:136,:);

% compute the values of yi:
N = data(:,8);  % speciment_7
A = data(:,1);

figure(1);
plot(A, N);
xlabel('A (Crack length)');
ylabel('N (Number of oscillations)');


%%
% Calculate dNdA Numerically:
dNdA = gradient(N,A); 

figure(2);
plot(A, dNdA);
xlabel('A (Crack length)');
ylabel('dN/dA');
legend('Numerical derivation');

%%
% regress Paris formula dN/dA = beta*A^alpha
% The Paris formula can be reexpressed as a linear equation: 
% ln(dN/dA) = alpha*ln(A)+ln(beta) <-> y = c1*x_1 + c_0*1  
% c_0 = ln(beta) --> beta = exp(c_0)
% c_1 = alpha

% Matlab
y = log(dNdA);
x_1 = log(A);
n = length(x_1);
X=[ones(n,1),x_1];
coeffs = regress(y,X);
c_0=coeffs(1);
c_1=coeffs(2);
y_regress =X*coeffs;

figure(3);
plot(x_1,y);
hold on;
correlation_y = corrcoef([x_1,y]);
title(['Correlation coefficient: ', num2str(correlation_y(1,end),'%0.2f')]);

plot(x_1,y_regress);
legend_regress = ['regression (c_0:',num2str(c_0,'%0.4f'), 'c_1:' num2str(c_1,'%0.4f'),')'];
legend('Numerical derivation',legend_regress);
xlabel('ln(A)');
ylabel('ln(dN/dA)');

beta = exp(c_0);
alpha = c_1;
dNdA_regress = beta.*A.^alpha;

figure(4);
plot(A, dNdA);
hold on;
xlabel('A (Crack length)');
ylabel('dN/dA');

plot(A, dNdA_regress);
legend_regress = ['regression (\alpha:',num2str(alpha,'%0.4f'), '\beta:' num2str(beta,'%0.0f'),')'];
legend('Numerical derivation',legend_regress);

N_regress = cumtrapz(A,dNdA_regress);
figure(5);
plot(A, N);
hold on;
plot(A, N_regress);
xlabel('A (Crack length)');
ylabel('N (Number of oscillations)');
legend('Measurement','Numerical integration of Paris Law','Location','NorthWest');

% Manual Least square fit:
error = y - (c_0 + c_1*x_1);

% First criteria: mean(error) = 0
% Second criteria: min(sum(error^2))
% (See solution in separate calculation)
S_y = sum(y);
S_x = sum(x_1);
As = y-S_y/n;
Bs = S_x/n-x_1;
c_1_ = -sum(As.*Bs)/sum(Bs.^2);
c_0_ = (S_y-c_1_*S_x)/n;

%%
% ANOVA
% Hypotesis: alpha=0
% --> c_1=0
% --> y = c_0

V1=1;
V2=(n-2);

X0=[ones(n,1)];
coeffs0 = regress(y,X0);
y_0 =X0*coeffs0;

figure(16);
plot(x_1,y);
hold on;
plot(x_1,y_regress);
plot(x_1,y_0);
legend('Data','Present regression model','H_0');
xlabel('ln(A)');
ylabel('ln(dN/dA)');
title('ANOVA test H_0:\alpha=0');

SSRu = sum((y-y_regress).^2);
SSRr = sum((y-y_0).^2);
F_ratio = ((SSRr-SSRu)/V1) / (SSRu/(V2));

F=finv(0.95,V1,V2);

figure(15);
X_ =linspace(0,2*F,100);
Y_ = fcdf(X_,V1,V2);
plot(X_,Y_);
hold on;
plot(F,0.95,'ro');
plot([F,F],[0,0.95],'k.:');
text(F,0.02,[' F=',num2str(F)]);
plot([0,F],[0.95,0.95],'k.:');
text(0,0.93,' P=0.95');

title(['F cumulative distribution function (DOF: ',num2str(V1),', ',num2str(V2), ')']);
xlabel('x');
ylabel('P(F<x)');

if F_ratio > F
    disp('alpha=0 Rejected!!!');
else
    disp('alpha=0 Accepted');
end

%%
% Regression diagnostics (check the residuals)

epsilon = y - y_regress;

% Check the normal distribution
figure(6)
normplot(epsilon);

% Check the constant variance
figure(7)
plot(y,epsilon, 'bo');
hold on;
xlabel('y=ln(dN/dA)');
ylabel('Residual error \epsilon');
correlation_y_error = corrcoef([y,epsilon]);
title(['Correlation coefficient: ', num2str(correlation_y_error(1,end),'%0.2f')]);
X2=[ones(n,1),y];
coeffs2 = regress(epsilon,X2);
epsilon_regress_y=X2*coeffs2;
plot(y,epsilon_regress_y, 'k:');
legend('\epsilon','error trend');


% Check the correlation between the residuals and independent variables
figure(8)
plot(x_1,epsilon,'r*');
hold on;
xlabel('x_1=ln(A)');
ylabel('Residual error \epsilon');
correlation_x_error = corrcoef([x_1,epsilon]);
title(['Correlation coefficient: ', num2str(correlation_x_error(1,end),'%0.2f')]);
X3=[ones(n,1),x_1];
coeffs3 = regress(epsilon,X3);
epsilon_regress_x_1=X2*coeffs3;
plot(x_1,epsilon_regress_x_1, 'k:');
legend('\epsilon','error trend');

n_parts = 2;
part_length = floor(length(epsilon)/n_parts);
figure(9);
hold on;
legends = {}
for nn=1:(n_parts)
    index_1 = 1+(nn-1)*part_length;
    index_2 = nn*part_length;
    x_= x_1(index_1:index_2);
    y_= epsilon(index_1:index_2);
    plot(x_, y_,'.');    
    legends{end+1} = ['var(\epsilon)=',num2str(var(y_))];
end
legend(legends);
xlabel('x_1=ln(A)');
ylabel('Residual error \epsilon');

%%
% Second part
%Compute the confidence intervals for   example   the   95%   confidence   intervals)  of   theestimated parameters ^α, ^β

[coeffs,coeffs_intervalls_95] = regress(y,X);
c_0_95_intervall = coeffs_intervalls_95(1,:);
c_1_95_intervall = coeffs_intervalls_95(2,:);
beta_95_intervall = exp(c_0_95_intervall);
alpha_95_intervall = c_1_95_intervall;

invall_alpha = 0.05;

%% c_0
c_0s = y-(c_1*x_1)
MU = mean(c_0s);
SIGMA = std(c_0s);
c_0_min = norminv(invall_alpha/2,MU,SIGMA);
c_0_max = norminv(1-invall_alpha/2,MU,SIGMA);
c_0_95_intervall2 = [c_0_min,c_0_max];

figure(10);
x_c0 = linspace(MU-3*SIGMA,MU+3*SIGMA,100);
Y = normcdf(x_c0,MU,SIGMA);

plot(x_c0,Y,'b-');
hold on;
plot([c_0_min,c_0_min],[0,normcdf(c_0_min,MU,SIGMA)],'k.:');
text(c_0_min,0.02,num2str(c_0_min));

plot([c_0_max,c_0_max],[0,normcdf(c_0_max,MU,SIGMA)],'k.:');
text(c_0_max,0.02,num2str(c_0_max));

plot([min(x_c0),c_0_min],invall_alpha/2*[1,1],'k.:');
text(min(x_c0),invall_alpha/2,num2str(invall_alpha/2));

plot([min(x_c0),c_0_max],(1-invall_alpha/2)*[1,1],'k.:');
text(min(x_c0),(1-invall_alpha/2),num2str((1-invall_alpha/2)));
title('Cumulative Normal distribution c_0')
xlabel('X');
ylabel('P(c_0<X)');

%%
%c_1
c_1s = (y-(c_0))./x_1;
MU = mean(c_1s);
SIGMA = std(c_1s);
c_1_min = norminv(invall_alpha/2,MU,SIGMA);
c_1_max = norminv(1-invall_alpha/2,MU,SIGMA);
c_1_95_intervall2 = [c_1_min,c_1_max];

figure(11);
x_c1 = linspace(MU-3*SIGMA,MU+3*SIGMA,100);
Y = normcdf(x_c1,MU,SIGMA);

plot(x_c1,Y,'b-');
hold on;
plot([c_1_min,c_1_min],[0,normcdf(c_1_min,MU,SIGMA)],'k.:');
text(c_1_min,0.02,num2str(c_1_min));

plot([c_1_max,c_1_max],[0,normcdf(c_1_max,MU,SIGMA)],'k.:');
text(c_1_max,0.02,num2str(c_1_max));

plot([min(x_c1),c_1_min],invall_alpha/2*[1,1],'k.:');
text(min(x_c1),invall_alpha/2,num2str(invall_alpha/2));

plot([min(x_c1),c_1_max],(1-invall_alpha/2)*[1,1],'k.:');
text(min(x_c1),(1-invall_alpha/2),num2str((1-invall_alpha/2)));
title('Cumulative Normal distribution c_1')
xlabel('X');
ylabel('P(c_1<X)');

%%
A = data(:,1);
alphas = [];
betas = [];
speciment_numbers = [1:6,8:68];
indexes = speciment_numbers+1;
for i=1:length(speciment_numbers)

    N_ = data(:,indexes(i));
    [alpha_, beta_] = regression_virkler(N_,A);
    alphas(i)=alpha_;
    betas(i)=beta_;
end

figure(12);
hist(alphas,10);
hold on;
plot(alpha_95_intervall,[0,0],'ro', 'MarkerSize',12,'MarkerFaceColor','r');
xlabel('\alpha');
ylabel('Number of values');
title('Histogram of \alpha for all the speciments');

figure(13);
hist(betas,10);
hold on;
plot(beta_95_intervall,[0,0],'ro', 'MarkerSize',12, 'MarkerFaceColor','r');

%%
% 99.99% Confidence interval.

[coeffs,coeffs_intervalls_9999] = regress(y,X,0.001);
c_0_9999_intervall = coeffs_intervalls_9999(1,:);
c_1_9999_intervall = coeffs_intervalls_9999(2,:);

beta_9999_intervall = exp(c_0_9999_intervall);
alpha_9999_intervall = c_1_9999_intervall;

figure(12);
plot(alpha_9999_intervall,[0,0],'gs', 'MarkerSize',12, 'MarkerFaceColor','g');
legend('All speciments','No 7 95% intervall', 'No 7 99.99% intervall');

figure(13);
plot(beta_9999_intervall,[0,0],'gs', 'MarkerSize',12, 'MarkerFaceColor','g');
legend('All speciments','No 7 95% intervall', 'No 7 99.99% intervall');

figure(14);
for i=1:length(speciment_numbers)
    N_ = data(:,indexes(i));
    p1 = plot(A,N_,'b-','LineWidth',2);
    hold on;
    p1.Color(4) = 0.15;
end;
xlabel('A (Crack length)');
ylabel('N (Number of oscillations)');
title('All the rest speciments');

%%
% Save figures:
exportgraphics(figure(1),'data.pdf');
exportgraphics(figure(2),'dNdA.pdf');
exportgraphics(figure(3),'ln_dNdA.pdf');
exportgraphics(figure(4),'dNdA_regression.pdf');
exportgraphics(figure(5),'N_regression.pdf');
exportgraphics(figure(6),'normplot.pdf');
exportgraphics(figure(7),'error_vs_y.pdf');
exportgraphics(figure(8),'error_vs_x.pdf');
exportgraphics(figure(9),'variance.pdf');
exportgraphics(figure(10),'distribution_c0.pdf');
exportgraphics(figure(11),'distribution_c1.pdf');
exportgraphics(figure(12),'all_alpha.pdf');
exportgraphics(figure(13),'all_beta.pdf');
exportgraphics(figure(14),'all_data.pdf');
exportgraphics(figure(15),'F-statistics.pdf');
exportgraphics(figure(16),'H0-regression.pdf');






