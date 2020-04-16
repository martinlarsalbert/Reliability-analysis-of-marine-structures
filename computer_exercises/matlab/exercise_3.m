close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo
%%
load('../Computer Exercise Data/testdata.mat')
figure();
plot(F,X,'bx')
