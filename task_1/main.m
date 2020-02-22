%% Параметры отрисовки
drawSet_n = 200; % число точек для отрисовки одного множества 
len = 6; % длина опорной гиперплоскости
%% Параметры алгоритма
psi_0_n = 50; % размер сетки для psi_0
n_check = 128; % количество направлений при проверке попадания в терминальное множество
h_t = 0.001; % шаг сетки по времени
T_max = 3; % максимально возможное t1
B_det_eps = 0.00001; % эпсилон для определителя В 


%% Описание входных данных
A = [1 0; 1 1] / 10;
B = [10 0; 0 10];
f = [1; 1];
a = 2;
b = 1;
x0 = [-3 7]';
r0 = 0.5;
x1 = [9 7 11 10;
      3 4 5 7];
eps = 1;

while det(B) < B_det_eps
    B = B + (rand(2) - 0.5) * 10 * sqrt(B_det_eps);
end

%% Вычисление траекторий, удовлетворяющих ПМП
clc;
phi = linspace(0, 2 * pi, psi_0_n + 1);
phi(psi_0_n + 1) = [];
psi_0 = [cos(phi); sin(phi)];

[X, U, Psi, opt_ind, t1, t] = calc_traj(psi_0, psi_0_n, A, B, f, h_t, T_max, a, b, x0, r0, x1, eps, n_check);

%% Проверка условия трансверсальности
x_1 = X(:, size(X, 2), opt_ind);
psi1 = Psi(:, size(Psi, 2), opt_ind);
psi1 = psi1 / norm(psi1);
l1 = find_dir(x_1, @(l) rho_X1(l, x1, eps), n_check);
angle = acos(-l1' * psi1) * 180 / pi;
disp("Transversality inaccuracy (angle in degrees) =");
disp(angle)

%% Уточнение вычислений
phi_step = phi(2) - phi(1);
phi = linspace(phi(opt_ind) - phi_step * psi_0_n / 10, phi(opt_ind) + phi_step * psi_0_n / 10, psi_0_n);
psi_0 = [cos(phi); sin(phi)];

[X, U, Psi, opt_ind, t1, t] = calc_traj(psi_0, psi_0_n, A, B, f, h_t, T_max, a, b, x0, r0, x1, eps, n_check);

%% График (x1, x2) - только ОУ
figure('Name', 'Plot (x1, x2)');
drawSet(@(l) rho_X0(l, x0, r0), drawSet_n, 'y');
hold on;
drawSet(@(l) rho_X1(l, x1, eps), drawSet_n, 'g');
xlabel('x_1');
ylabel('x_2');
axis equal;

if opt_ind
    plot(X(1, :, opt_ind), X(2, :, opt_ind), 'r', 'LineWidth', 2);
end

x_0 = X(:, 1, opt_ind);
psi0 = psi_0(:, opt_ind);
l = x_0 - x0;
l = l / norm(l);
tau = [-l(2); l(1)];
plot(x_0(1) + linspace(-len/2, len/2, drawSet_n) * tau(1), x_0(2) + linspace(-len/2, len/2, drawSet_n) * tau(2), 'black', 'LineWidth', 2);
hold on;
plot(x_0(1) + linspace(0, 1, drawSet_n) * psi0(1), x_0(2) + linspace(0, 1, drawSet_n) * psi0(2), 'black', 'LineWidth', 2);
hold on;

x_1 = X(:, size(X, 2), opt_ind);
psi1 = Psi(:, size(Psi, 2), opt_ind);
psi1 = psi1 / norm(psi1);
l1 = find_dir(x_1, @(l) rho_X1(l, x1, eps), n_check);
tau = [-l1(2); l1(1)];
plot(x_1(1) + linspace(-len/2, len/2, drawSet_n) * tau(1), x_1(2) + linspace(-len/2, len/2, drawSet_n) * tau(2), 'black', 'LineWidth', 2);
hold on;
plot(x_1(1) + linspace(0, 1, drawSet_n) * (-psi1(1)), x_1(2) + linspace(0, 1, drawSet_n) * (-psi1(2)), 'black', 'LineWidth', 2);
hold on;

%% График (x1, x2) - Все траектории
figure('Name', 'Plot (x1, x2)');
drawSet(@(l) rho_X0(l, x0, r0), drawSet_n, 'y');
hold on;
drawSet(@(l) rho_X1(l, x1, eps), drawSet_n, 'g');
xlabel('x_1');
ylabel('x_2');
axis equal;

for i = 1:psi_0_n
    plot(X(1, :, i), X(2, :, i), 'b');
    hold on;
end

if opt_ind
    plot(X(1, :, opt_ind), X(2, :, opt_ind), 'r', 'LineWidth', 2);
end

x_0 = X(:, 1, opt_ind);
psi0 = psi_0(:, opt_ind);
l = x_0 - x0;
l = l / norm(l);
tau = [-l(2); l(1)];
plot(x_0(1) + linspace(-len/2, len/2, drawSet_n) * tau(1), x_0(2) + linspace(-len/2, len/2, drawSet_n) * tau(2), 'black', 'LineWidth', 2);
hold on;
plot(x_0(1) + linspace(0, 1, drawSet_n) * psi0(1), x_0(2) + linspace(0, 1, drawSet_n) * psi0(2), 'black', 'LineWidth', 2);
hold on;

x_1 = X(:, size(X, 2), opt_ind);
psi1 = Psi(:, size(Psi, 2), opt_ind);
psi1 = psi1 / norm(psi1);
l1 = find_dir(x_1, @(l) rho_X1(l, x1, eps), n_check);
tau = [-l1(2); l1(1)];
plot(x_1(1) + linspace(-len/2, len/2, drawSet_n) * tau(1), x_1(2) + linspace(-len/2, len/2, drawSet_n) * tau(2), 'black', 'LineWidth', 2);
hold on;
plot(x_1(1) + linspace(0, 1, drawSet_n) * (-psi1(1)), x_1(2) + linspace(0, 1, drawSet_n) * (-psi1(2)), 'black', 'LineWidth', 2);
hold on;
%% График (t, x1)
figure('Name', 'Plot (t, x1)');
xlabel('t');
ylabel('x1');
axis equal;

for i = 1:psi_0_n
    plot(t, X(1, :, i), 'b');
    hold on;
end

if opt_ind
    plot(t, X(1, :, opt_ind), 'r', 'LineWidth', 2);
end

%% График (t, x2)
figure('Name', 'Plot (t, x2)');
xlabel('t');
ylabel('x2');
axis equal;

for i = 1:psi_0_n
    plot(t, X(2, :, i), 'b');
    hold on;
end

if opt_ind
    plot(t, X(2, :, opt_ind), 'r', 'LineWidth', 2);
end

%% График (u1, u2)
figure('Name', 'Plot (u1, u2)');
drawSet(@(l) rho_P(l, a, b), drawSet_n, 'b');
hold on;
xlabel('u_1');
ylabel('u_2');
plot(U(1, :, opt_ind), U(2, :, opt_ind), 'or', 'MarkerSize', 10);
hold on;
axis equal;

%% График (t, u1)
figure('Name', 'Plot (t, u1)');
xlabel('t');
ylabel('u1');
axis equal;

for i = 1:psi_0_n
    plot(t, U(1, :, i), 'b');
    hold on;
end

if opt_ind
    plot(t, U(1, :, opt_ind), 'r', 'LineWidth', 2);
end

%% График (t, u2)
figure('Name', 'Plot (t, u2)');
xlabel('t');
ylabel('u2');
axis equal;

for i = 1:psi_0_n
    plot(t, U(2, :, i), 'b');
    hold on;
end

if opt_ind
    plot(t, U(2, :, opt_ind), 'r', 'LineWidth', 2);
end

%% График (psi1, psi2)
figure('Name', 'Plot (psi1, psi2)');
xlabel('psi_1');
ylabel('psi_2');

for i = 1:psi_0_n
    plot(Psi(1, :, i), Psi(2, :, i), 'b');
    hold on;
end

if opt_ind
    plot(Psi(1, :, opt_ind), Psi(2, :, opt_ind), 'r', 'LineWidth', 2);
end

%% График (t, psi1)
figure('Name', 'Plot (t, psi1)');
xlabel('t');
ylabel('psi1');
axis equal;

for i = 1:psi_0_n
    plot(t, Psi(1, :, i), 'b');
    hold on;
end

if opt_ind
    plot(t, Psi(1, :, opt_ind), 'r', 'LineWidth', 2);
end

%% График (t, psi2)
figure('Name', 'Plot (t, psi2)');
xlabel('t');
ylabel('psi2');
axis equal;

for i = 1:psi_0_n
    plot(t, Psi(2, :, i), 'b');
    hold on;
end

if opt_ind
    plot(t, Psi(2, :, opt_ind), 'r', 'LineWidth', 2);
end