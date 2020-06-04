close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo
%%
W=[5.3, 8.1, 5.5, 4.0, 6.1, 9.1, 10.8, 7.0, 5.2, 8.0, 9.3, 8.4, 6.5, 4.5, 3.8, 2.9, 3.7, 4.1, 3.5, 5.7, 4.3, 3.0, 3.5, 4.2, 4.0, 2.8];
index = 1:length(W);

mask = W(1:end-2)<W(2:end-1) & W(2:end-1)>W(3:end);
index_mask = find(mask)+1;
Wp = W(index_mask);
POT_mask = Wp > 5;
index_POT = index_mask(find(POT_mask));
Wpot=W(index_POT);

figure(1);
plot(index,W,'.-');
hold on;
plot(index_mask,Wp,'x');
plot(index_POT,Wpot,'o');




