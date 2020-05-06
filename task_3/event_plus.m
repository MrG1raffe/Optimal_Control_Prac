function [value, isterminal, direction] = event_plus(~, x)
%     value = [x(4); max(x(2), 5*x(1)^5 + x(1) * sin(x(1)^3) - 1)];
%     isterminal = [1; 0];
%     direction = [-1; 0];
    value = x(4);
    isterminal = 1;
    direction = -1;
end