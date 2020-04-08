function [value, isterminal, direction] = event(~, y, psi0, umax, M, dir)
    value = [y(3) - 4 * psi0 * y(2) * umax^3; y(3); y(2) - M];
    isterminal = [1; 1; 1];
    direction = [dir; -1; -1];
end