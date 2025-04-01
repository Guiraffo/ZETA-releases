demoID = 'Zonotope sampling: uniform, Gaussian, grid';
util.demoinit;

Z = Zonotope( zeros(2,1), [ 1, 0,  1;
                            1, 2, -1]);
           
nof_samples = 1000;
grid_divs = 30;
SampleUnif = sample(Z,nof_samples);
SampleGauss = sample(Z,nof_samples,'Gaussian');
SampleGrid = sample(Z,grid_divs,'grid');  


figure
plot(Z,'blue',0); 
hold on;
plot(SampleUnif(1,:),SampleUnif(2,:),'b.');
xlabel('x_1');
ylabel('x_2');
title('Uniform distribution: 1000 samples');

figure
plot(Z,'blue',0); 
hold on;
plot(SampleGauss(1,:),SampleGauss(2,:),'r.');
xlabel('x_1');
ylabel('x_2');
title('Truncated Gaussian distribution: 1000 samples');
           
figure
plot(Z,'blue',0); 
hold on;
plot(SampleGrid(1,:),SampleGrid(2,:),'k.');
xlabel('x_1');
ylabel('x_2');
title('Grid sampling: 20 divisions');
