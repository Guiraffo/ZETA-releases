demoID = 'Parallelotope enclosing a constrained zonotope';
util.demoinit;

Z = CZonotope([0.0353; -0.8925; -0.0427], [ 0.4994,    0.9011,    0.4358,   -0.3923,   -0.6777,    0.3984,    0.9341;
                                           -0.3173,    0.2396,    0.5523,    0.4305,   -0.6426,    0.6147,   -0.7389;
                                           -0.6745,    0.4580,    0.8033,   -0.5168,   -0.3678,   -0.7603,    0.1995],...
                                          [-0.7243,   -0.2214,   -0.2121,   -0.3954,   -0.6819,   -0.7319,   -0.6745;
                                           -0.1650,   -0.1881,    0.0173,    0.3305,   -0.1544,    0.3794,    0.6618;
                                            0.8288,   -0.4967,   -0.9112,    0.9685,   -0.7604,    0.6913,   -0.1271],...
                                          [-0.4708;    0.2866;    0.7428]);
                                      
Zhull = intervalhull(Z);

Zparallelred = partopebound(Z,'red');
Zparallelresc = partopebound(Z,'rescale');


%% Volumes (requires MPT)
Zpoly = Polyhedron(Z);
Zhullpoly = Polyhedron(Zhull);
Zparallelredpoly = Polyhedron(Zparallelred);
Zparallelrescpoly = Polyhedron(Zparallelresc);
disp(strcat(['Volume of Z: ',num2str(volume(Zpoly))]))
disp(strcat(['Volume of Zhull: ',num2str(volume(Zhullpoly))]))
disp(strcat(['Volume of Zparallelred: ',num2str(volume(Zparallelredpoly))]))
disp(strcat(['Volume of Zparallelresc: ',num2str(volume(Zparallelrescpoly))]))


%% Figures

figure
hold on;
plot(Z,'blue',1);
plot(Zhull,'yellow',0);
plot(Zparallelred,'yellow',0);
plot(Zparallelresc,'red',0.3);
xlabel('x1');
ylabel('x2');
zlabel('x3');



