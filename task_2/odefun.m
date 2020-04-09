function f = odefun(~, y, psi30, psi0, g, umax, mode)
    if (mode == 1)
        u = umax;
    end
    if (mode == 2)
        u = nthroot(y(3) / (4 * psi0 * y(2)), 3);
    end
    if (mode == 3)
        u = 0;
    end
    f = zeros(4, 1);
    f(1) = -g + u * y(1) / y(2);
    f(2) = -u;
    f(3) = -psi30 * y(1) - y(4) * g - y(3) * u / y(2);
    f(4) = -psi30 - y(4) * u / y(2);
end