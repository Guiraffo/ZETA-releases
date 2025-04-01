demoID = 'Propagating a zonotope through a nonlinear function';
util.demoinit;

%% Inputs

Z = Interval([0.1; 0.2], [0.9;0.8]);
fname = 'Model1';

%% Settings

OPTIONS_IANat = util.propagopts('IANat');
OPTIONS_IAMV  = util.propagopts('IAMV');

%% Propagate using different approaches

Z_IANat = propagate(Z,fname,OPTIONS_IANat);
Z_IAMV = propagate(Z,fname,OPTIONS_IAMV);


%% Samples

Z0_samp = sample(Z,20,'grid');
nof_samp = size(Z0_samp,2);
Z_samp = zeros(size(Z,1),nof_samp);
for j=1:nof_samp
    Z_samp(:,j) = feval(str2func(fname),Z0_samp(:,j));
end


%% Figures

figure
plot(Z,'yellow',1);
xlabel('x_1');
ylabel('x_2');
title('Z to propagate')
grid off;
box on;

figure
hold on;
plot(Z_IAMV,'blue',0.7);
plot(Z_IANat,'red',0.7);
plot(intersect(Z_IAMV,Z_IANat),'magenta',1);
plot(Z_samp(1,:),Z_samp(2,:),'k.');
xlabel('x_1');
ylabel('x_2');
title('f(Z) using IANat and IAMV')
legend('IANat','IAMV','IANat \cap IAMV')
grid off;
box on;

