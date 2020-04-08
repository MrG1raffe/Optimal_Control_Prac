function [x_opt, u_opt, tau_switch_opt, H] = calc_optimal_traj_1(T, M, m0, umax, l, g, eps, ~)
    % Calculates an optimal control, optimal trajectory 
    % and switch time for problem 1
    
    time_N = 1000;
    tau_N = 1000;
    t = linspace(0, T, time_N);
    if (eps >= l)
        error("Incorrect input. l must be greater than eps.");
    end
    if (umax < m0 * g / l)
        error("Incorrect input (umax < m0 * g / l)");
    end
    if (T <= 0 || M <= 0 || m0 <= 0 || umax <= 0 || l <= 0 || g <= 0 || eps <= 0)
        error("Incorrect input. All values must be positive.");
    end
    if (m0 <= M)
        error("Incorrect input. m0 must be greater than M.");
    end
    
    tau_f = (m0 - M) / umax;
    if (tau_f  < T && x1_umax(T, M, m0, umax, l, g) < l - eps)
        disp("The problem has no solution.");
        x1 = x1_umax(t, M, m0, umax, l, g);
        plot(t, x1 - l, 'b', t, (-eps) * ones(size(t)), 'r--', t, (eps) * ones(size(t)), 'r--');
        xlabel('t');
        ylabel('v');
        return;
    end
    tau_s = linspace(min(tau_f, T), 0, tau_N);
    x1_T = x1_umax(tau_s, M, m0, umax, l, g) - g * (T - tau_s);
    idx_opt = find(x1_T <= l + eps, 1, 'first');
    tau_switch_opt = tau_s(idx_opt);
    t = linspace(0, T, time_N);
    x1 = x1_umax(t, M, m0, umax, l, g, tau_switch_opt);
    x2 = max(m0 - umax * t, m0 - umax * tau_switch_opt);
    x_opt = [x1; x2];
    u_opt = umax * (t <= tau_switch_opt);
    x3 = height(t, x_opt(1, :), l);
    H = x3(numel(x3));
    disp("Switch time:");
    disp(tau_switch_opt);
    disp("H(T) = ");
    disp(H);
    if (nargin == 8)
          plot_res(t, x1, x2, x3, u_opt, l, eps);
    end
end