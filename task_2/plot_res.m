function plot_res(t, x1, x2, x3, u, l, eps, M, m0, umax)
    T = t(numel(t));
    
    subplot(2, 2, 1);
    plot(t, x1 - l, 'b', t, (-eps) * ones(size(t)), 'r--', t, (eps) * ones(size(t)), 'r--');
    xlabel('t');
    ylabel('v');
    xlim([0 T]);
    ylim([- 2 * eps, max(x1 - l) + eps]);
    %axis equal;

    subplot(2, 2, 2);
    plot(t, x2, 'b', t, M * ones(size(t)), 'r--',  t, m0 * ones(size(t)), 'r--');
    xlabel('t');
    ylabel('m');
    xlim([0 T]);
    %axis equal;

    subplot(2, 2, 4);
    plot(t, u, 'b', t, zeros(size(t)), 'r--', t, umax * ones(size(t)), 'r--');
    xlabel('t');
    ylabel('u');
    xlim([0 T]);

    subplot(2, 2, 3);
    plot(t, x3, 'b', t, zeros(size(t)), 'r--');
    xlabel('t');
    ylabel('H');
    xlim([0 T]);
    %axis equal;   
end