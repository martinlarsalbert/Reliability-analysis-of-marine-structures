function [alpha, beta] = regression_virkler(N,A)

% Calculate dNdA Numerically:
dNdA = gradient(N,A); 

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