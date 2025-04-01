demoID = 'Rescaling a constrained zonotope';
util.demoinit;


%% Input constrained zonotope

Z = CZonotope( zeros(2,1), [ 1, 0,  1;
                             1, 2, -1],...
               [-2, 1, -1], 2);
        
%% Rescale using linear programming and interval arithmetic           

ZrescaledLP = rescale(Z,'LP');
ZrescaledIA = rescale(Z,'IA');


%% Figures

figure
hold on;
plot(Z,'blue',0.8);
plot(Zonotope(Z.c,Z.G),'yellow',0);
plot(Zonotope(ZrescaledLP.c,ZrescaledLP.G),'yellow',0.2);
title('Linear programming')
xlabel('x_1');
ylabel('x_2');

figure
hold on;
plot(Z,'blue',0.8);
plot(Zonotope(Z.c,Z.G),'yellow',0);
plot(Zonotope(ZrescaledIA.c,ZrescaledIA.G),'yellow',0.2);
title('Interval arithmetic')
xlabel('x_1');
ylabel('x_2');
