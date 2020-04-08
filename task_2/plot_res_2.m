function plot_res_2(t, x1, x2, x3, u, l, H, umax, M)
    T = t(numel(t));
    
    subplot(2, 2, 1);
    plot(t, x1 - l, 'b');
    xlabel('t');
    ylabel('v');
    %axis equal;

    subplot(2, 2, 2);
    y = x2;
    plot(t, y, 'b', t, M * ones(size(t)), 'r--');
    xlabel('t');
    ylabel('m');
    %axis equal;

    subplot(2, 2, 4);
    plot(t, u, 'b', t, umax * ones(size(t)), 'r--');
    xlabel('t');
    ylabel('u');
    xlim([0 T]);
    ylim([0 umax]);

    subplot(2, 2, 3);
    plot(t, x3, 'b', t, H * ones(size(t)), 'r--');
    xlabel('t');
    ylabel('H');
    %axis equal;   
end