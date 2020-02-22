function [val, point] = rho_P(l, a, b)
    if a < 1
        b = b / a;
        a = 1 / a;
    end
    if abs(l(1)) > abs(l(2)) && abs(l(2)) * sqrt(a * (a + 1)) < sqrt(l(1)^2 + a * l(2)^2)
        val = sqrt(b * (l(1)^2 + a*l(2)^2) / a);
        point = [l(1) * sqrt(b / a), l(2) * sqrt(a * b)]' / sqrt(l(1)^2 + a * l(2)^2);
    else
        if abs(l(2)) > abs(l(1)) && abs(l(1)) * sqrt(a * (a + 1)) < sqrt(a * l(1)^2 + l(2)^2)
            val = sqrt(b * (a * l(1)^2 + l(2)^2) / a);
            point = [l(1) * sqrt(b * a), l(2) * sqrt(b / a)]' / sqrt(a * l(1)^2 + l(2)^2);
        else
            val = sqrt(b / (a + 1)) * (abs(l(1)) + abs(l(2)));
            point = sqrt(b / (a + 1)) * [sign(l(1)), sign(l(2))]; 
        end
    end
end