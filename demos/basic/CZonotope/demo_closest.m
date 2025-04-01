demoID = 'Closest point in a constrained zonotope or in its interval hull to a desired point in space';
util.demoinit;

%% Inputs

Z = CZonotope( zeros(2,1), [ -2, 0,  1;
                              1, 2, -1],...
               [-2, 1, -1], 2);
         
%% Operations           
           
Zhull = intervalhull(Z);           

z1 = Z.c;
z2 = [3.5;3];
z3 = [-1.5;1];
z4 = [1;3];
z5 = [2;-2];

h1_Z = closest(Z,z1,'set');
h2_Z = closest(Z,z2,'set');
h3_Z = closest(Z,z3,'set');
h4_Z = closest(Z,z4,'set');
h5_Z = closest(Z,z5,'set');

h1_hull = closest(Z,z1,'hull');
h2_hull = closest(Z,z2,'hull');
h3_hull = closest(Z,z3,'hull');
h4_hull = closest(Z,z4,'hull');
h5_hull = closest(Z,z5,'hull');

%% Figures

figure;
hold on;
plot(Z,'blue',0.3);
plot(z1(1),z1(2),'rx')
plot(h1_Z(1),h1_Z(2),'r+')
plot(z2(1),z2(2),'gx')
plot(h2_Z(1),h2_Z(2),'g+')
plot(z3(1),z3(2),'mx')
plot(h3_Z(1),h3_Z(2),'m+')
plot(z4(1),z4(2),'bx')
plot(h4_Z(1),h4_Z(2),'b+')
plot(z5(1),z5(2),'kx')
plot(h5_Z(1),h5_Z(2),'k+')
axis([-2 4 -3 4])
box on;
grid on;

figure;
hold on;
plot(Z,'blue',0.3);
plot(Zhull,'blue',0);
plot(z1(1),z1(2),'rx')
plot(h1_hull(1),h1_hull(2),'r+')
plot(z2(1),z2(2),'gx')
plot(h2_hull(1),h2_hull(2),'g+')
plot(z3(1),z3(2),'mx')
plot(h3_hull(1),h3_hull(2),'m+')
plot(z4(1),z4(2),'bx')
plot(h4_hull(1),h4_hull(2),'b+')
plot(z5(1),z5(2),'kx')
plot(h5_hull(1),h5_hull(2),'k+')
axis([-2 4 -3 4])
box on;
grid on;
