close all;
clear all;
clc;
addpath('../wafo_old')
initwafo
%%
% Random variables:
% W : section modulus
% M : Load

mu_w = 39;  % [m3]
sigma_w = 3.2;
mu_m = 6.2e9; % [Nm]
sigma_m = 1.1e9 % [Nm]
sigma_a = 160e6; % Max stress [Pa]
gamma = 1.5; % Safety factor

% m = r - s  % limit state function
% where:
% r : resistance
% s : load
% resistance is M=sigma*w (and times gamma in this case)
% r = sigma_a*W*gamma
% s = M

mu_r = sigma_a*mu_w*gamma
mu_s = mu_m
sigma_r = sigma_a*sigma_w*gamma
sigma_s = sigma_m

beta_c = (mu_r-mu_s) / sqrt(sigma_r^2 + sigma_s^2)

