demoID = 'Convex hull of a collection of constrained zonotopes';
util.demoinit;

%% Input sets
   
X = CZonotope([5; 5],  [1,  1,  -1;
                        0,  1,   1], [1, -1, 1], -1);

W = CZonotope(zeros(2,1)+5, [-1, 1,  0;
                              2, 0, -1], [1, -1, 1], -1);
                        
Y = CZonotope(-ones(2,1), [-1,  1,  1;
                           -1, -1,  2], [1, -1, 1], -1);                        
                        
%% Convex hull of two and three CZs                        

convXW = convexhull(X,W);
convXWY = convexhull(X,W,Y);

%% Figures

figure
hold on
plot(convXW,'magenta',0.8,'Yalmip');
plot(X,'green',0.8);
plot(W,'blue',0.8);
axis tight;
xlabel('x_1')
ylabel('x_2')
title('Convex hull of two CZs');
legend('Convex hull','X','W');
grid off;
box on;

figure
hold on
plot(convXWY,'magenta',0.8,'Yalmip');
plot(X,'green',0.8);
plot(W,'blue',0.8);
plot(Y,'yellow',0.8);
axis tight
xlabel('x_1')
ylabel('x_2')
title('Convex hull of three CZs');
legend('Convex hull','X','W','Y');
grid off;
box on;

