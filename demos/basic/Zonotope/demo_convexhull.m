demoID = 'Convex hull of a collection of zonotopes';
util.demoinit;

%% Input sets
   
X = Zonotope([5; 5],  [1,  1,  -1;
                       0,  1,   1]);

W = Zonotope(zeros(2,1)+5, [-1, 1,  0;
                             2, 0, -1]);
                        
Y = Zonotope(-ones(2,1), [-1,  1,  1;
                          -1, -1,  2]);                        
                        
%% Convex hull of two and three CZs                        

convXW = convexhull(X,W);
convXWY = convexhull(X,W,Y);

%% Figures

figure
hold on
plot(convXW,'magenta',0.8);
plot(X,'green',0.8);
plot(W,'blue',0.8);
axis tight;
xlabel('x_1')
ylabel('x_2')
title('Convex hull of two zonotopes');
legend('Convex hull','X','W');

figure
hold on
plot(convXWY,'magenta',0.8);
plot(X,'green',0.8);
plot(W,'blue',0.8);
plot(Y,'yellow',0.8);
axis tight
xlabel('x_1')
ylabel('x_2')
title('Convex hull of three zonotopes');
legend('Convex hull','X','W','Y');

