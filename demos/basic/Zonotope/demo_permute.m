demoID = 'Permute the dimensions of a zonotope';
util.demoinit;



Z = Zonotope( [0; -1], [ 1, 0,  1;
                          1, 2, -1]);
           

Zperm = permute(Z,[2,1],[]);

figure;
hold on
plot(Z,'green',0.3);
plot(Zperm,'yellow',0.3);
legend('Original','Permuted')
xlabel('x_1');
ylabel('x_2');