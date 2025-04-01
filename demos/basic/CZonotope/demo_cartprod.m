demoID = 'Cartesian product of two constrained zonotopes';
util.demoinit;


%% Inputs

Z = CZonotope( zeros(2,1), [ 1, 0,  1;
                             1, 2, -1],...
               [-2, 1, -1], 2);           
Y = CZonotope( 1, 0.5);

%% Operations
ZxY = [Z;Y];


%% Command window display and figures

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
box on;
grid off;
figure;
plot(Y,'blue',0.3);
title('Y');
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
box on;
grid off;
figure;
plot(ZxY,'magenta',0.3);
title('Z \times Y');
xlabel('x_1'); ylabel('x_2'); zlabel('x_3');
box on;
grid off;



