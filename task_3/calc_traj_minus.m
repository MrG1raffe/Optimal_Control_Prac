function [end_of_traj, sws] = calc_traj_minus(alpha, tspan, x0, sw_cur, num_of_sw)
    if num_of_sw <= 3
        sw_cur(num_of_sw, 1) = x0(1);
        sw_cur(num_of_sw, 2) = x0(2);
    end
    opt = odeset('Events', @event_minus);
    [~, x, te, ye, ~] = ode45(@(t, x) odefun_minus(t, x, alpha), tspan, x0, opt);
    %plot(x(:, 1), x(:,2), 'y');
    hold on;
    if isempty(te)
        end_of_traj = [x(end, 1), x(end, 2)];
        sws = sw_cur;
    else
        [end_of_traj, sws] = calc_traj_plus(alpha, [te, tspan(end)], ye, sw_cur, num_of_sw + 1);
    end
end