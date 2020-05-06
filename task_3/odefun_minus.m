function f = odefun_minus(~, x, alpha)
    f = zeros(4, 1);
    f(1) = x(2);
    f(2) = -alpha - 5 * x(1).^5 - x(1) .* sin(x(1).^3) - x(2);
    f(3) = x(4) .* (25 * x(1).^4 + sin(x(1).^3) + 3 * (x(1).^3) .* cos(x(1).^3));
    f(4) = x(4) - x(3);
end