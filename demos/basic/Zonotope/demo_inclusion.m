demoID = 'Zonotope inclusion';
util.demoinit;

p = [1;1];

M = [Interval(2,3), Interval(1,1.5), 1, Interval(-3,-2.8); Interval(-1,-0.5), Interval(0,1), Interval(3,4), Interval(5,6)];

Z = Zonotope.inclusion(p,M);

Zinf = Zonotope(p,inf(M));
Zmid = Zonotope(p,mid(M));
Zsup = Zonotope(p,sup(M));


figure;
hold on;
plot(Z,'magenta',1);
plot(Zsup,'green',1);
plot(Zmid,'red',1);
plot(Zinf,'blue',1);
legend('inclusion','Zsup','Zmid','Zinf')