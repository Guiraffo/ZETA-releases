demoID = 'Parallelotope enclosing of a zonotope';
util.demoinit;


G = [-0.2812, -0.1968,  0.4235;
      0.5186, -0.2063, -0.2267];
 
Z = Zonotope(zeros(2,1), G);
   

Zhull = intervalhull(Z);
Zpartope = partopebound(Z);

figure
hold on
plot(Zhull,'cyan',0.1);
plot(Zpartope,'yellow',0.5);
plot(Z,'magenta',0.9);
xlabel('x_1');
ylabel('x_2');
legend('Interval hull','Parallelotope enclosure','Z');