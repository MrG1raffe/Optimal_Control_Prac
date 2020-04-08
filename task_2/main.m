%% Входные данные
T = 2;
M = 1;
m0 = 3;
umax = 1;
H = 0.5;
l = 1;
g = 0.2;
eps = 0.01;

%% Задача 1
clc;
[~, ~, ~, ~] = calc_optimal_traj_1(T, M, m0, umax, l, g, eps, 'plot');

%% Задача 2
clc;
[x, u, t_switch, J_min] = calc_optimal_traj_2(T, M, m0, umax, l, g, H, 'plot');


%% Пример режимов 2-1-3
T = 2;
M = 1;
m0 = 6;
umax = 4;
H = 2.3;
l = 1;
g = 0.5;
eps = 0.01;