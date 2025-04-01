demoID = 'Division of two intervals using polyhedral relaxations.';
util.demoinit;

%% Inputs and operations


xinterval = Interval(-1,1);
yinterval = Interval(1,2);
[xsample,ysample] = meshgrid(inf(xinterval):0.05:sup(xinterval),inf(yinterval):0.05:sup(yinterval));
plots.color = 'yellow';

Polyrelax.clear;

X = Polyrelax(xinterval);
Y = Polyrelax(yinterval);
Z = X/Y;
zsample = xsample./ysample;

Hrep = Polyrelax.Hrep;
thisZ = Polyrelax.Z;
[thisZpoly.H,thisZpoly.k] = hrep(thisZ);

% Division requires intersection with original interval for reduced
% conservatism
[HrepNthisZ.H,HrepNthisZ.k] = util.intersectpoly(Hrep.H,Hrep.k,[],[],thisZpoly.H,thisZpoly.k,[],[]);

%% Figures

figure;
surf(xsample,ysample,zsample);
hold on;
util.plotpoly(Hrep.H,Hrep.k,plots.color,0.1);
title('z = x/y')
xlabel('x')
ylabel('y')
zlabel('z')
box on;
grid off;

figure;
surf(xsample,ysample,zsample);
hold on;
util.plotpoly(HrepNthisZ.H,HrepNthisZ.k,plots.color,0.1);
title('z = x/y, \cap Z')
xlabel('x')
ylabel('y')
zlabel('z')
box on;
grid off;