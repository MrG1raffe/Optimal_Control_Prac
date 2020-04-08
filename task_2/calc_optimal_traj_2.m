function [x_opt, u_opt, tau_s_opt, J_min] = calc_optimal_traj_2(T, M, m0, umax, l, g, H, ~)
    % Calculates an optimal control, optimal trajectory, 
    % moments of switches and minimal value of J for problem 2     
    
    time_N = 1000;
    H_eps = 0.05 * H;
    psi1_eps = 0.5;
    psi_N = 40;
    
    x_opt = [];
    u_opt = [];
    tau_s_opt = [];
    J_min = [];
    
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
    if (abs(x3(end) - H) < H_eps)
        if (tau_f < T)
            tau_s_opt = tau_f;
        end
        u_opt = umax * (t <= tau_f);
        J_min = trapz(t, u_opt.^4);
        x2 = max(m0 - umax * t, M);
        x_opt = [x1; x2];
        if nargin == 8
            plot_res_2(t, x1, x2, x3, u_opt, l, H, umax, M);
        end
        return;
    end
    
    J_min = umax^4 * min(tau_f, T);
    modes_opt = [];
    
    psi10 = linspace(-1, 1, 2 * psi_N);
    psi30 = linspace(-1, 1, 2 * psi_N);
    psi0 = linspace(0, 1, psi_N);
    psi0(1) = [];
    psi_min_dif = 100;
    H_min_dif = 100;
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
                 disp("Calculating trajectory");
                while (te < tspan(2))
                    [t, y, te, ye, ie] = ode45(@(t, y) odefun(t, y, psi30(j), psi0(k), g, umax, mode), tspan, y0, options);
                    if ~isempty(te)
                        te = te(end);
                        ye = ye(size(ye, 1), :);
                        ie = ie(end);
                    end
                    y_res = [y_res; y];
                    t_res = [t_res; t];
                    y0 = ye;
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
                if (abs(psi1(end)) < psi_min_dif)
                    psi_min_dif = abs(psi1(end));
                end
                if (abs(psi1(end)) < psi1_eps)
                    x3 = height(t_res, x1, l);
                    if abs(x3(end) - H) < H_min_dif
                        H_min_dif = abs(x3(end) - H);
                    end
                    if (abs(x3(end) - H) < H_eps)
                        
                        u = (F > 4 .* x2 * psi0(k) * umax^3) * umax + (F <= 4 * x2 * psi0(k) * umax^3) .* (F > 0) .* nthroot(F ./ (4 * x2 * psi0(k)), 3) + 0;
                        u = u .* (t_res <= t_turn_off);
                        J = trapz(t_res, u.^4);
                        if (J < J_min)
                            psi30_opt = psi30(j);
                            t_opt = t_res;
                            u_opt = u;
                            x_opt = [x1; x2; x3];
                            F_opt = F;
                            tau_s_opt = switches;
                            J_min = J;
                            psi0_opt = psi0(k);
                            modes_opt = modes;
                        end
                    end
                end
            end
        end
    end
    psi30_opt
    if isempty(modes_opt)
        H_min_dif
        psi_min_dif
        disp("Optimal control wasn't found.");
    else
        disp("Optimal control was found.");
        disp("J_min = ");
        disp(J_min);
        disp("Modes:");
        disp(modes_opt);
        disp("Switch times:");
        disp(tau_s_opt);
        figure
        plot(t_opt, F_opt, t_opt, F_opt - 4 * psi0_opt * x_opt(2,:) * umax^3);
        legend('F', 'K');
        if nargin == 8
            figure
            plot_res_2(t_opt, x_opt(1, :), x_opt(2, :), x_opt(3, :), u_opt, l, H, umax, M);
        end
    end
end
    