demoID = 'Zonotope intersection with strips';
util.demoinit;

G = [0.2812,  0.1968,  0.4235;
     0.0186, -0.2063, -0.2267];

Z = Zonotope(zeros(2,1), G);

S = Strip([1; -1], 0, 0.1);

Zinter1 = intersection(Z,S,'seg');
Zinter2 = intersection(Z,S,'Bravo');


figure
hold on
plot(Z,'magenta',1);
plot(S,'black',0);
axis0 = axis;
plot(Zinter1,'green',1);
plot(Zinter2,'yellow',1);
axis(axis0);
xlabel('x_1');
ylabel('x_2');
hlegend = legend('Z','S','Z \cap S (seg)','Z \cap S (Bravo)');