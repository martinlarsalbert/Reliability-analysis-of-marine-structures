function [a_new, beta_new, k] = hasofer(a,beta)

gamma = 35*a(1)*a(2)*beta + 350*a(1) + 350*a(2) -300*a(3);
beta_new = -2000/gamma;

df = [35*a(2)*beta+350;
      35*a(1)*beta+350;
      -300];

k=sqrt(sum(df.^2)); 
a_new = -df/k;