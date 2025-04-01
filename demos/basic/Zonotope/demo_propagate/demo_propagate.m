demoID = 'Propagating a zonotope through a nonlinear function';
util.demoinit;

%% Inputs

Z = Zonotope([0.5;0.5],[0.2, 0.1,  0.1;
                          0, 0.2, -0.1]);
fname = 'Model1';

%% Settings

OPTIONS_ZMV = util.propagopts('ZMV');
OPTIONS_ZFO = util.propagopts('ZFO');
OPTIONS_ZDC = util.propagopts('ZDC');

%% Propagate using Mean Value Extension and First-Order Taylor Extension

Z_MV = propagate(Z,fname,OPTIONS_ZMV);
Z_FO = propagate(Z,fname,OPTIONS_ZFO);
Z_DC = propagate(Z,fname,OPTIONS_ZDC);

%% Samples

Z0_samp = sample(Z,20,'grid');
nof_samp = size(Z0_samp,2);
Z_samp = zeros(Z.dim,nof_samp);
for j=1:nof_samp
    Z_samp(:,j) = feval(str2func(fname),Z0_samp(:,j));
end


%% Figures

figure
plot(Z,'yellow',1);
xlabel('x_1');
ylabel('x_2');
title('Z to propagate')


figure
hold on;
plot(Z_MV,'blue',1);
plot(Z_DC,'cyan',1);
plot(Z_FO,'red',1);
plot(Z_samp(1,:),Z_samp(2,:),'k.');
xlabel('x_1');
ylabel('x_2');
title('f(Z) using MV, FO, DC')
legend('MV','DC','FO')

