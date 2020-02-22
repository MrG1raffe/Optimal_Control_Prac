function res = is_in_terminal(point, rho_set, n)
    phi = linspace(0, 2 * pi, n + 1);
    phi(n + 1) = [];
    phi = phi(randperm(n));
    l = [cos(phi); sin(phi)];
    res = 1;
    for i = 1:n
        [val, ~] = rho_set(l(:, i));
        if l(:, i)' * point > val
            res = 0;
            break
        end
    end
end