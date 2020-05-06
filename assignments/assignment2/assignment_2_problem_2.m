close all;
clear all;
clc;
addpath('../../wafo_old')
initwafo

%%
data_course_1 = load('COURSE1.mat');  % mean values and covarianve matrix for course 1
data_course_2 = load('COURSE2.mat');  % mean values and covarianve matrix for course 1

N=100;
result_course_1 = simulate(data_course_1,N);

result_course_2 = simulate(data_course_2,N);

%%
% Save figures:
%exportgraphics(figure(1),'data.pdf');







