a = 3;
x1 = linspace(-1/sqrt(a), 1/sqrt(a), 1000);
y1 = sqrt(1 - a*x1.^2);
plot(x1, y1, 'b', x1, -y1, 'b')
axis equal;
hold on
x2 = linspace(-1, 1, 1000);
y2 = sqrt(1 - x2.^2) / sqrt(a);
plot(x2, y2, 'b', x2, -y2, 'b')
hold on
c = sqrt(1/(a+1));
xp = linspace(-c, c, 1000);
yp = sqrt(1 - xp.^2) / sqrt(a);
plot(xp,yp,'r',...
    'LineWidth',3)
plot(xp,-yp,'r',...
    'LineWidth',3)
plot(yp,xp,'r',...
    'LineWidth',3)
plot(-yp,xp,'r',...
    'LineWidth',3)
axis([-1.2 1.2 -1.2 1.2]);