function [val, point] = rho_X1(l, points, eps)
    [mx, ind] = max(l' * points);
    val = mx + eps * norm(l);
    point = eps * l / norm(l) + points(:, ind);
end