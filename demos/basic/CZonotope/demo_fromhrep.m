demoID = 'Conversion of a convex polytope in halfspace representation into CG-rep';
util.demoinit;

%% Inputs

H1 = [ 1,  2;
      -1, -2;
       1,  1;
      -1, -1];
k1 = [1; 1; 1; 1]; 

H2 = [ 1,  2;
      -1, -2;
       1,  1;
      -1, -1;
       1, -1;
      -1,  1];
k2 = [1; 1; 1; 1; 2; 2]; 



H3 = [ 1, -2,  1;
      -1,  2, -1;
       1,  1, -1;
      -1, -1,  1;
       1,  1,  1;
      -1, -1, -1];  
k3 = [1; 1; 1; 1; 1; 1]; 

%% Operations

Z = CZonotope.fromhrep(H1,k1);
Z2 = CZonotope.fromhrep(H2,k2);
Z3 = CZonotope.fromhrep(H3,k3);


%% Figures

figure
util.plotpoly(H1,k1,'black',0.3);
title('P_1 in H-rep');
xlabel('x1');
ylabel('x2');

figure;           
plot(Z,'yellow',0.3);
title('P_1 in CG-rep');
xlabel('x1');
ylabel('x2');

figure
util.plotpoly(H2,k2,'black',0.3);
title('P_2 in H-rep');
xlabel('x1');
ylabel('x2');

figure;           
plot(Z2,'yellow',0.3);
title('P_2 in CG-rep');
xlabel('x1');
ylabel('x2');

figure
util.plotpoly(H3,k3,'black',0.3);
title('P_3 in H-rep');
xlabel('x1');
ylabel('x2');
zlabel('x3');

figure;           
plot(Z3,'yellow',0.3);
title('P_3 in CG-rep');
xlabel('x1');
ylabel('x2');
zlabel('x3');
       