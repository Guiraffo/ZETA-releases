demoID = 'Constrained zonotope sampling: uniform, grid, pseudo-uniform';
util.demoinit;

Z = CZonotope( zeros(2,1), [ 1, 0,  1;
                             1, 2, -1],...
               [-2, 1, -1], 2);
           
nof_samples = 1000;
grid_divs = 30;
SampleUnif = sample(Z,nof_samples);
SampleGrid = sample(Z,grid_divs,'grid');  
SamplePUnif = sample(Z,nof_samples);


figure
plot(Z,'blue',0); 
hold on;
plot(SampleUnif(1,:),SampleUnif(2,:),'b.');
xlabel('x_1');
ylabel('x_2');
title('Uniform distribution: 1000 samples');
grid off;
box on;
           
figure
plot(Z,'blue',0); 
hold on;
plot(SampleGrid(1,:),SampleGrid(2,:),'k.');
xlabel('x_1');
ylabel('x_2');
title('Grid sampling: 20 divisions');
grid off;
box on;

figure
plot(Z,'blue',0); 
hold on;
plot(SamplePUnif(1,:),SamplePUnif(2,:),'r.');
xlabel('x_1');
ylabel('x_2');
title('Pseudo-uniform distribution: 1000 samples');
grid off;
box on;