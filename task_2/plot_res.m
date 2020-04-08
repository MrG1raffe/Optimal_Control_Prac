function plot_res(t, x1, x2, x3, u, l, eps)
    T = t(numel(t));
    
    subplot(2, 2, 1);
    plot(t, x1 - l, 'b', t, (-eps) * ones(size(t)), 'r--', t, (eps) * ones(size(t)), 'r--');
    xlabel('t');
    ylabel('v');
    ylim([- 2 * eps, max(x1 - l) + eps]);
    %axis equal;

    subplot(2, 2, 2);
    y = x2;
    plot(t, y, 'b');
    xlabel('t');
    ylabel('m');
    %axis equal;

    subplot(2, 2, 4);
    plot(t, u, 'b');
    xlabel('t');
    ylabel('u');
    xlim([0 T]);

    subplot(2, 2, 3);
    plot(t, x3, 'b');
    xlabel('t');
    ylabel('H');
    %axis equal;   
end