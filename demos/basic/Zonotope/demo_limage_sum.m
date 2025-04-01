demoID = 'Linear image and Minkowsky sum of zonotopes';
util.demoinit;


A = [ 0.2, -0.1;
      0.5,  0.2];
 
B = [0.5,    1;
       0,  0.5];
   

X0 = Zonotope([15; 15],  [3,  0;
                           0,  5]);

W = Zonotope(zeros(2,1), 1*eye(2));

nof_iterations = 2;



X = cell(nof_iterations+1,1);
% Xinterval = cell(nof_iterations+1,1);

X{1} = X0;
% Xinterval{1} = midrad(X0.c,diag(X0.G));

for k=1:nof_iterations
    X{k+1} = A*X{k} + B*W;
%     Xinterval{k+1} = A*Xinterval{k} + B*midrad(W.c,diag(W.G));
end


figure;
hold on;
plot(X0,'green',0.4);
for k=2:nof_iterations+1
    plot(X{k},'yellow',0.4);
%     plot_box(Xinterval{k},'gray',0);
end
xlabel('x_1')
ylabel('x_2')

