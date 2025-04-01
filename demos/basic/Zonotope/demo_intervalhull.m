demoID = 'Zonotope interval hull';
util.demoinit;

G = [0.2812,  0.1968,  0.4235;
     0.0186, -0.2063, -0.2267];

Z = Zonotope(zeros(2,1), G);
Zhull = intervalhull(Z);

figure
hold on
plot(Z,'magenta',0.7);
plot(Zhull,'yellow',0.1);
xlabel('x_1');
ylabel('x_2');
