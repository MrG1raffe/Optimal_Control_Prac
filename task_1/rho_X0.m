function [val, point] = rho_X0(l, x0, r0)
    point = x0 + r0 * l / norm(l);
    val = l' * point;
end