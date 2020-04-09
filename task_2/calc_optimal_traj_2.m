function [t_opt,x_opt, u_opt, tau_s_opt, J_min] = calc_optimal_traj_2(T, M, m0, umax, l, g, H, psi_N, H_eps, ~)
    % Calculates an optimal control, optimal trajectory, 
    % moments of switches and minimal value of J for problem 2     
    
    time_N = 1000;
    %H_eps = 0.01 * H;
    psi1_eps = 0.5;
    %psi_N = 200;
    
    t_opt = [];
    x_opt = [];
    u_opt = [];
    tau_s_opt = [];
    J_min = [];
    psi_opt = [];
    
    if (umax < m0 * g / l)
        error("Incorrect input (umax < m0 * g / l)");
    end
    if (T <= 0 || M <= 0 || m0 <= 0 || umax <= 0 || l <= 0 || g <= 0 || H <= 0)
        error("Incorrect input. All values must be positive.");
    end
    if (m0 <= M)
        error("Incorrect input. m0 must be greater than M.");
    end
    
    t = linspace(0, T, time_N);
    x1 = x1_umax(t, M, m0, umax, l, g);
    x3 = height(t, x1, l);
    if (x3(end) < H)
        disp("The problem has no solution.");
        plot(t, x3, 'b', t, (H) * ones(size(t)), 'r--');
        xlabel('t');
        ylabel('H');
        return;
    end
    
    tau_f = (m0 - M) / umax;
    tau_s_opt = [];
    if (abs(x3(end) - H) < 0.001)
        if (tau_f < T)
            tau_s_opt = tau_f;
        end
        u_opt = umax * (t <= tau_f);
        J_min = trapz(t, u_opt.^4);
        x2 = max(m0 - umax * t, M);
        x_opt = [x1; x2];
        if nargin == 10
            plot_res_2(t, x1, x2, x3, u_opt, l, H, umax, M, m0);
        end
        return;
    end
    
    J_min = umax^4 * T;
    modes_opt = [];
    
    psi10 = linspace(-1, 1, 2 * psi_N);
    psi30 = linspace(-1, 1, 2 * psi_N);
    psi0 = linspace(0, 1, psi_N);
    psi0(1) = [];
    H_min_dif = inf;
    for i=1:(2 * psi_N)
        for j = 1:(2 * psi_N)
            if (psi10(i) * psi30(j) < 0 || abs(psi30(j)) * T > abs(psi10(i)) || psi10(i)^2 + psi30(j)^2 > 1)
                continue;
            end
            for k = 1:(psi_N - 1)
                if (psi10(i)^2 + psi30(j)^2 + psi0(k)^2 > 1)
                    continue;
                end
                psi20 = -sqrt(1 - psi10(i)^2 - psi30(j)^2 - psi0(k)^2);
                F0 = psi10(i) * l - psi20 * m0;            
                if (F0 < 4 * psi0(k) * m0^4 * g^3 / l^3)
                    continue;
                end
                if (F0 > 4 * m0 * psi0(k) * umax^3)
                    mode = 1;
                else
                    mode = 2;
                end
                if (psi30(j) > 0)
                    dir = -1;
                else
                    dir = 1;
                end
                options = odeset('Events', @(t, x) event(t, x, psi0(k), umax, M, dir));
                y0 = [l; m0; F0; psi10(i)];
                tspan = [0 T];
                te = tspan(1);
                y_res = [];
                t_res = [];
                switches = [];
                modes = mode;
                t_turn_off = T;
                str = ['Calculating trajectory: psi0 = ', num2str(psi0(k)), ' psi10 = ', num2str(psi10(i)), ' psi20 = ', num2str(psi20), ' psi30 = ', num2str(psi30(j))];
                disp(str);
                while (~isempty(te))
                    [t, y, te, ye, ie] = ode45(@(t, y) odefun(t, y, psi30(j), psi0(k), g, umax, mode), tspan, y0, options);
                    if ~isempty(te)
                        te = te(end);
                        ye = ye(size(ye, 1), :);
                        ie = ie(end);
                    end
                    y_res = [y_res; y];
                    t_res = [t_res; t];
                    y0 = ye;
                    if (t(end) == T)
                        break;
                    end
                    if (~isempty(ie)) && (mode ~= 3)
                        if (ie == 1)
                            if (dir == -1)
                                mode = 2;
                            else
                                mode = 1;
                            end
                            
                        end
                        if ((ie == 2) || (ie == 3))
                            if (ie == 3)
                                t_turn_off= te;
                            end
                            mode = 3;
                        end
                        if (~isempty(modes) && mode ~= modes(end))
                            modes = [modes, mode];
                        end
                    end
                    if ~isempty(te)
                        tspan(1) = te;
                        if (isempty(switches) || (~isempty(switches) && (te - switches(end)) > T * 0.001 && t_turn_off == T))
                            switches = [switches, te];
                        end
                    end
                end
                
                t_res = t_res';
                x1 = y_res(:, 1)';
                x2 = y_res(:, 2)';
                F = y_res(:, 3)';
                psi1 = y_res(:, 4)';
                if (abs(psi1(end)) < psi1_eps)
                    x3 = height(t_res, x1, l);
                    if 100 * abs(x3(end) - H) / H < H_min_dif
                        H_min_dif = 100 * abs(x3(end) - H);
                    end
                    if (abs(x3(end) - H) < H_eps)
                        u = (F > 4 .* x2 * psi0(k) * umax^3) * umax + (F <= 4 * x2 * psi0(k) * umax^3) .* (F > 0) .* nthroot(F ./ (4 * x2 * psi0(k)), 3) + 0;
                        u = u .* (t_res <= t_turn_off);
                        J = trapz(t_res, u.^4);
                        if (J < J_min)
                            t_opt = t_res;
                            u_opt = u;
                            x_opt = [x1; x2; x3];
                            tau_s_opt = switches;
                            J_min = J;
                            psi_opt = [psi0(k), psi10(i), psi20, psi30(j)]; 
                            modes_opt = modes;
                        end
                    end
                end
            end
        end
    end
    if isempty(modes_opt)
        disp("Optimal control wasn't found.");
        disp("Minimal height inaccuracy, %:");
        disp(H_min_dif);
    else
        disp("Optimal control was found.");
        disp("Optimal parameters: psi0 = ");
        disp(psi_opt);
        disp("J_min = ");
        disp(J_min);
        disp("Modes:");
        disp(modes_opt);
        disp("Switch time(s):");
        disp(tau_s_opt);
%         figure
%         plot(t_opt, F_opt);
%         legend('F);
        if nargin == 10
            figure
            plot_res_2(t_opt, x_opt(1, :), x_opt(2, :), x_opt(3, :), u_opt, l, H, umax, M, m0);
        end
    end
end
    