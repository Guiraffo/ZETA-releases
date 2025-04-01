demoID = 'Zonotope support strips across given directions';
util.demoinit;

rng('default');
rng(1250);

Z  = Zonotope(zeros(2,1), [ -0.1200, -0.1000,  0.5000;
                             0.2000, -0.1000, -0.1000]);

r1 = [5; 5];
r2 = [2;-4];

S1 = support_strip(Z,r1);
S2 = support_strip(Z,r2);
     


figure
hold on;
plot(Z,'green',1);
plot(S1,'red',0)
plot(S2,'blue',0)
axis0 = axis;
plot(Z,'green',1);
axis(axis0*2);
xlabel('x_1');
ylabel('x_2');