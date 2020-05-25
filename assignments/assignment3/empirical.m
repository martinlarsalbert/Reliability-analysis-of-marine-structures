function [x,y] = empirical(x_s) 
N=length(x_s); 
x=sort(x_s);
y=(1:N)/N;