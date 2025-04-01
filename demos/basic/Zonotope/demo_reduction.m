demoID = 'Zonotope generator reduction';
util.demoinit;


G = [0, -0.1200, 0, -0.0100, 0,  0.0500, -1.5000, -0.7500, 0, -0.4500, 0, -1.0380, 0,       0;
     0,  0.0200, 0, -0.1000, 0, -0.1100,  1.5000,  0.0000, 0,  0.4500, 0,  2.0760, 0,  2.2788];

Z = Zonotope(zeros(2,1), G);
ngmax = 3;

Zreduc1 = reduction(Z,ngmax,'Combastel');
Zreduc2 = reduction(Z,ngmax,'CombastelW',diag([0.5,1]));
Zreduc3 = reduction(Z,ngmax,'Girard');
Zreduc4 = reduction(Z,ngmax,'Chisci');

figure
hold on;
plot(Zreduc4,'red',1);
plot(Z,'blue',1)
xlabel('x_1');
ylabel('x_2');
legend('Reduced Z: Chisci','Z');

figure
hold on
plot(Zreduc2,'magenta',0.8);
plot(Zreduc1,'green',0.8);
plot(Zreduc3,'yellow',0.8);
plot(Zreduc4,'red',0.8);
plot(Z,'blue',1)
xlabel('x_1');
ylabel('x_2');
legend('CombastelW','Combastel','Girard','Chisci','Z');