demoID = 'Pontryagin difference between a convex polytope and a zonotope';
util.demoinit;

Z = Zonotope( zeros(2,1), [ 1, 0,  1;
                             1, 2, -1]);

B = Interval([-2;-1],[0;2]);

Bzon = Zonotope(mid(B),diag(rad(B)));

[H_z,k_z] = hrep(Z);

[H_p,k_p] = pontdiff(Bzon,H_z,k_z);

figure;
plot(Z,'green',0.8);
hold on
plot(B,'red',0.7)
plot(Polyhedron(H_p,k_p),'Color','yellow','Alpha',0.7);
xlabel('x_1');
ylabel('x_2');
legend('Z','B','Z-B')



