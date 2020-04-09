%% Входные данные
T = 5;
M = 1;
m0 = 2;
umax = 2;
H = 20;
l = 1;
g = 0.2;
eps = 0.01;

%% Задача 1
clc;
[~, ~, ~, ~] = calc_optimal_traj_1(T, M, m0, umax, l, g, eps, 'plot');

%% Задача 2
clc;
psi_N = 60;
H_eps = H * 0.1;
[t, x, u, t_switch, J_min] = calc_optimal_traj_2(T, M, m0, umax, l, g, H, psi_N, H_eps, 'plot');
