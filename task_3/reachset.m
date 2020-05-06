function [X, Y, Sx, Sy] = reachset(alpha, T, t_size)
    max_sw = 3;
    X = zeros(2 * t_size, 1);
    Y = zeros(2 * t_size, 1);
    Sx = zeros(max_sw, 2 * t_size) * nan;
    Sy = Sx;

    tspan = [0 T];
    
    opt = odeset('Events', @event_x2_plus);
    [~, ~, te, ~, ~] = ode45(@(t, x) odefun_plus(t, x, alpha), tspan, [0 0 0 0]', opt);
    if ~isempty(te)
        tau_plus = te;
    else
        tau_plus = T;
    end
    
    opt = odeset('Events', @event_x2_minus);
    [~, ~, te, ~, ~] = ode45(@(t, x) odefun_minus(t, x, alpha), tspan, [0 0 0 0]', opt);
    if ~isempty(te)
        tau_minus = te;
    else
        tau_minus = T;
    end


    
    t_plus = linspace(0, min(T, tau_plus), t_size + 1);
    [~, x_plus] = ode45(@(t, x) odefun_plus(t, x, alpha), t_plus, [0 0 0 0]'); 
    t_plus(1) = [];
    for i = 1:t_size
        %i
        tspan = [t_plus(i), T];
        if t_plus(i) == T                
            X(t_size) = x_plus(end, 1);
            Y(t_size) = x_plus(end, 2);
            break;
        end
        x0 = [x_plus(i + 1, 1); x_plus(i + 1, 2); 1; 0];
        [tr_end, sws] = calc_traj_minus(alpha, tspan, x0, nan*zeros(max_sw, 2), 1);
        Sx(:, i) = sws(:, 1);
        Sy(:, i) = sws(:, 2);
        X(i) = tr_end(1);
        Y(i) = tr_end(2);
    end

    
    t_minus = linspace(0, min(T, tau_minus), t_size + 1);
    [~, x_minus] = ode45(@(t, x) odefun_minus(t, x, alpha), t_minus, [0 0 0 0]'); 
    t_minus(1) = [];
    for i = 1:t_size
        tspan = [t_minus(i), T];
        if t_minus(i) == T                
            X(end) = x_minus(end, 1);
            Y(end) = x_minus(end, 2);
            break;
        end
        x0 = [x_minus(i + 1, 1); x_minus(i + 1, 2); -1; 0];
        [tr_end, sws] = calc_traj_plus(alpha, tspan, x0, nan*zeros(max_sw, 2), 1);
        Sx(:, i + t_size) = sws(:, 1);
        Sy(:, i + t_size) = sws(:, 2);
        X(i + t_size) = tr_end(1);
        Y(i + t_size) = tr_end(2);
    end
    [X, Y] = normalize(X, Y);
    X = [X, X(1)];
    Y = [Y, Y(1)];
end