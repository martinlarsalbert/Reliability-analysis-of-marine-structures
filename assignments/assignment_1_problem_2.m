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

y = log(dNdA);
x_1 = log(A);
X=[ones(length(x_1),1),x_1];
coeffs = regress(y,X);
c_0=coeffs(1);
c_1=coeffs(2);

beta = exp(c_0);
alpha = c_1;
dNdA_regress = beta.*A.^alpha;
plot(A, dNdA_regress);
legend_regress = ['regress \alpha:',num2str(alpha,'%0.4f'), '\beta:' num2str(beta,'%0.0f')];
legend('Numerical derivation',legend_regress);