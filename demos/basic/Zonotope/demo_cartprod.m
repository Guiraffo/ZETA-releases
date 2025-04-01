demoID = 'Cartesian product of two zonotopes';
util.demoinit;


Z = Zonotope( zeros(2,1), [ 1, 0, 1; 1, 2, -1]);
Y = Zonotope( 1, 0.5);
ZxY = [Z,Y];

disp('Z.c = '); disp(Z.c);
disp('Z.G = '); disp(Z.G);
disp('Y.c = '); disp(Y.c);
disp('Y.G = '); disp(Y.G);
disp('ZxY.c = '); disp(ZxY.c);
disp('ZxY.G = '); disp(ZxY.G);

figure;
plot(Z,'red',0.3);
title('Z');
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
figure;
plot(Y,'blue',0.3);
title('Y');
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
figure;
plot(ZxY,'magenta',0.3);
title('Z \times Y');
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');



