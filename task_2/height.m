function x3 = height(t, x1, l)
    % returns vector of rocket's height
    x3 = cumtrapz(t, x1 - l);
end