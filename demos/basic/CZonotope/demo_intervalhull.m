demoID = 'Interval hull of a constrained zonotope';
util.demoinit;

Z = CZonotope( zeros(2,1), [ 1, 0,  1;
                             1, 2, -1],...
               [-2, 1, -1], 2);
           
      
Wc = [-0.5704; -0.4547; -0.2853];
WG = [ 0.2372,  0.1605,  0.2727, -0.4427;
       0.2827, -0.3322,  0.0108,  0.3819;
       0.4280,  0.5464, -0.7632, -0.6770];
WA = [0.7536, -0.0318,  0.9153, -0.7388];
Wb = 0.2691;

W = CZonotope(Wc, WG, WA, Wb);

Zbox = intervalhull(Z);
Wbox = intervalhull(W);

figure;
plot(Z,'blue',1);
hold on;
plot(Zbox,'yellow',0.1);
xlabel('x_1');
ylabel('x_2');

figure
plot(W,'blue',0.7);
hold on
plot(Wbox,'yellow',0.1)
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');