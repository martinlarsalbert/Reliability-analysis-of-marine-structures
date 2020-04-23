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
figure(2);
plot(A, dNdA_regress);
legend_regress = ['regression (\alpha:',num2str(alpha,'%0.4f'), '\beta:' num2str(beta,'%0.0f'),')'];
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

xest = y_regress;

SSR=sum((xest-mean(y)).^2);
SSE=sum((y-xest).^2);
SST=sum((y-mean(y)).^2);
% you could check if SST=SSR+SSE
if ((SST-(SSR+SSE))/(SSR+SSE)>0.01)
    warning('SST=SSR+SSE Failed');
end

MSR = SSR/1;
MSE = SSE/(n-2);

F_ratio=(MSR)/(MSE);
P=1-fcdf(F_ratio, 1,n-2);

if F_ratio > P
    disp('alpha=0 Rejected!!!');
else
    disp('alpha=0 Accepted');
end

%%
% Regression diagnostics (check the residuals)

epsilon = y - y_regress;

% Check the normal distribution
figure()
normplot(epsilon);

% Check the constant variance
figure()
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
figure()
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
figure();
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




