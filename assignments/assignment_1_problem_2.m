close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo

%%
load('Virkler.mat');  % Loading into Virk matrix

% Note that only the first 136 rows could be used for the analysis in thisp roject. 
data = Virk(1:136,:);

% compute the values of yi:
N = data(:,8);  % speciment_7
A = data(:,1);

figure();
plot(A, N);
xlabel('A (Crack length)');
ylabel('N (Number of oscillations)');

%%
% Calculate dNdA Numerically:
dNdA = gradient(N,A); 



figure();
plot(A, dNdA);
xlabel('A (Crack length)');
ylabel('dN/dA');
hold on;

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

beta = exp(c_0);
alpha = c_1;
dNdA_regress = beta.*A.^alpha;
plot(A, dNdA_regress);
legend_regress = ['regress \alpha:',num2str(alpha,'%0.4f'), '\beta:' num2str(beta,'%0.0f')];
legend('Numerical derivation',legend_regress);

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
y_estimation = c_0;
error_estimation = y - y_estimation;
X=error_estimation;
SSR=sum((error_estimation-mean(error_estimation)).^2);
SSE=sum((y_estimation-error_estimation).^2);
SST=sum((y_estimation-mean(y_estimation)).^2);
% you could check if SST=SSR+SSE
F_ratio=(SSR/1)/(SSE/7);
P=1-fcdf(x_1, 1,7);





