
%%
t_size = 100;
T = 100;
alpha = 0.000001;
[X, Y, Sx, Sy] = reachset(alpha, T, t_size);

for i = 1:3
    if mod(i, 2)
        plot(Sx(i,1:t_size), Sy(i,1:t_size), 'r', Sx(i,t_size + 1:2*t_size), Sy(i,t_size + 1:2*t_size), 'b');
        hold on;
    else
        plot(Sx(i,1:t_size), Sy(i,1:t_size), 'b', Sx(i,t_size + 1:2*t_size), Sy(i,t_size + 1:2*t_size), 'r');
        hold on;
    end
end

plot(X, Y, 'black');

%% Мультик
mov = reachsetdyn(0.8, 0, 8, 50, 'homevideo.avi');

%%
movie(mov, 1, 5);