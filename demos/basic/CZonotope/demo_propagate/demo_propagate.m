demoID = 'Propagating a constrained zonotope through a nonlinear function';
util.demoinit;

%% Inputs

Z = CZonotope([-1;1],[0.2, 0.4,  0.2;
                      0.2,   0, -0.2], [2, 2, 2], -3);
fname = 'Model1';

%% Settings

OPTIONS_CZMV = util.propagopts('CZMV');
OPTIONS_CZFO = util.propagopts('CZFO');
OPTIONS_CZDC = util.propagopts('CZDC');
OPTIONS_CZPR = util.propagopts('CZPR');

%% Propagate using Mean Value Extension and First-Order Taylor Extension

Z_MV = propagate(Z,fname,OPTIONS_CZMV);
Z_FO = propagate(Z,fname,OPTIONS_CZFO);
Z_DC = propagate(Z,fname,OPTIONS_CZDC);
Z_PR = propagate(Z,fname,OPTIONS_CZPR);

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
plot(Z_FO,'red',1);
plot(Z_MV,'blue',1);
plot(Z_DC,'cyan',1);
plot(Z_PR,'magenta',1,'fast');
plot(Z_samp(1,:),Z_samp(2,:),'k.');
xlabel('x_1');
ylabel('x_2');
title('f(Z) using MV, FO, DC, PR')
legend('FO','MV','DC','PR')
grid off;
box on;

