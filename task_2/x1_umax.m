function x1 = x1_umax(t, M, m0, umax, l, g, tau_s)
    % returns value x1(t) with control u(t) = umax with engine 
    % turn off at moment tau_s
    tau_f = (m0 - M) / umax;
    if (nargin == 6)
        tau_s = tau_f;
    end
   
    x1 = g * m0 / (2 * umax) - g * min(t, tau_s) / 2 + m0 * (g * m0 - 2 * l * umax) ./ (2 * umax * (-m0 + min(t, tau_s) * umax));
    x1 = x1 - g * (t - tau_s) .* (t > tau_s);
end