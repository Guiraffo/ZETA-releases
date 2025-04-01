demoID = 'Intersection of a constrained zonotope and a convex polytope';
util.demoinit;

%% 2D

Z = CZonotope( zeros(2,1), [ 1, 0,  1;
                             1, 2, -1],...
               [-2, 1, -1], 2);
           
[H_p,k_p] = hrep(Interval(-ones(2,1),zeros(2,1)));           


ZnP = intersection(Z,H_p,k_p);

figure
hold on
plot(Z,'blue',1);
util.plotpoly(H_p,k_p,'red',1);
plot(ZnP,'magenta',1);
xlabel('x_1')
ylabel('x_2')
legend('Z','P','Z \cap P');
title('Intersection of a constrained zonotope and a convex polytope')
