function mov = reachsetdyn(alpha, t1, t2, n, filename)
    t_size = 100;
    mov(1:n) = struct('cdata', [], 'colormap', []);
    if ~t1
        T = linspace(t1, t2, n+1);
        T(1) = [];
    else
        T = linspace(t1, t2, n);
    end
    for j = 1:n
        j
        [X, Y, Sx, Sy] = reachset(alpha, T(j), t_size);
        plot(X, Y, 'black');
        hold on;
        for i = 1:3
            if mod(i, 2)
                plot(Sx(i,1:t_size), Sy(i,1:t_size), 'r', Sx(i,t_size + 1:2*t_size), Sy(i,t_size + 1:2*t_size), 'b');
                hold on;
            else
                plot(Sx(i,1:t_size), Sy(i,1:t_size), 'b', Sx(i,t_size + 1:2*t_size), Sy(i,t_size + 1:2*t_size), 'r');
                hold on;
            end
        end
        axis([-1.1 1.1 -1.2 1.2])
        mov(j) = getframe;
        cla;
    end
    
    if nargin == 5
        v = VideoWriter(filename);
        open(v);
        writeVideo(v, mov);
        close(v);
    end
end