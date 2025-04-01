demoID = 'Polyhedral relaxations of nonlinear functions.';
util.demoinit;

%%

disp('Warning: This demo requires MPT.')

% Default values
funcflags.all = 0; % Show all examples
funcflags.inv = 0; % 1/x
funcflags.exp = 0; % exp(x)
funcflags.log = 0; % ln(x)
funcflags.pow = 0; % x^a
funcflags.sin = 0; % sin(x)
funcflags.cos = 0; % cos(x)
funcflags.tan = 0; % tan(x)
funcflags.sec = 0; % sec(x)
funcflags.sqr = 0; % sqrt(x)
funcflags.abs = 0; % abs(x)
funcflags.nrm = 0; % norm(x,2) and norm(x,1)
funcflags.pni = 0; % non-integer power

% Setup only the ones of interest
funcflags.all = 1;


j = 0;

plots.color = 'y';
plots.grid = 'OFF';

%% 1/(ln(x))
if(funcflags.all||funcflags.inv)

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(5,10));
tic; Z{j} = 1/X{j}; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = 1./xsample{j};
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
P6int = P{j} & Polyhedron(Polyrelax.Z);
plot(P6int.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = 1/x')
xlabel('x')
ylabel('z')
grid(plots.grid);


end
%% exp(x)
if(funcflags.all||funcflags.exp)
    
j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-1,3));
tic; Z{j} = exp(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = exp(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = exp(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-1,3));
tic; Z{j} = exp(-X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = exp(-xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = exp(-x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

end
%% log(x)
if(funcflags.all||funcflags.log)
    
j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(1,10));
tic; Z{j} = log(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = log(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = log(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(1,10));
tic; Z{j} = -log(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = -log(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = -log(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

end
%% 1/(ln(x))
if(funcflags.all||funcflags.log)

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(5,10));
tic; Z{j} = 1/(log(X{j})); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = 1./(log(xsample{j}));
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
P6int = P{j} & Polyhedron(Polyrelax.Z);
plot(P6int.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = 1/ln(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

end
%% x^a (a even)
if(funcflags.all||funcflags.pow)
    
j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-2,4));
tic; Z{j} = X{j}^2; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^2;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^2')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-5,0));
tic; Z{j} = X{j}^2; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^2;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^2')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(1,8));
tic; Z{j} = X{j}^2; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^2;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^2')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-8,1));
tic; Z{j} = X{j}^2; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^2;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^2')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-2,1));
tic; Z{j} = X{j}^6; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^6;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^6')
xlabel('x')
ylabel('z')
grid(plots.grid);

end
%% x^a (a odd)
if(funcflags.all||funcflags.pow)
    
j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(2,4));
tic; Z{j} = X{j}^3; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^3;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^3')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-4,-1));
tic; Z{j} = X{j}^3; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^3;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^3')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-4,5));
tic; Z{j} = X{j}^3; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^3;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^3')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-2,1));
tic; Z{j} = X{j}^3; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^3;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^3')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-5,10));
tic; Z{j} = X{j}^3; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^3;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^3')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-2,3));
tic; Z{j} = X{j}^5; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^5;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^5')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-2,3));
tic; Z{j} = X{j}^7; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^7;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^7')
xlabel('x')
ylabel('z')
grid(plots.grid);

end
%% sin(x)
if(funcflags.all||funcflags.sin)
    
j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(4*pi/2,5*pi/2));
tic; Z{j} = sin(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = sin(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sin(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(3.5*pi/2,4.7*pi/2));
tic; Z{j} = sin(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = sin(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sin(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(5.5*pi/2,6.5*pi/2));
tic; Z{j} = sin(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = sin(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sin(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(5*pi/2,6.9*pi/2));
tic; Z{j} = sin(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = sin(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sin(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(3.1*pi/2,6.9*pi/2));
tic; Z{j} = sin(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = sin(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sin(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(4.1*pi/2,6.9*pi/2));
tic; Z{j} = sin(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = sin(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sin(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(0,pi/2));
tic; Z{j} = sin(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = sin(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sin(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-pi,pi));
tic; Z{j} = sin(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = sin(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sin(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);


j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-2.1*pi-1.5*pi,1.9*pi/2));
tic; Z{j} = sin(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = sin(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sin(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

end
%% cos(x)
if(funcflags.all||funcflags.cos)
    
j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(4*pi/2,5*pi/2));
tic; Z{j} = cos(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = cos(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = cos(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(3.5*pi/2,4.7*pi/2));
tic; Z{j} = cos(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = cos(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = cos(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(5.5*pi/2,6.5*pi/2));
tic; Z{j} = cos(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = cos(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = cos(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(5*pi/2,6.9*pi/2));
tic; Z{j} = cos(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = cos(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = cos(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(3.1*pi/2,6.9*pi/2));
tic; Z{j} = cos(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = cos(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = cos(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);


j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(4.1*pi/2,6.9*pi/2));
tic; Z{j} = cos(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = cos(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = cos(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(0,pi/2));
tic; Z{j} = cos(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = cos(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = cos(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-pi,pi));
tic; Z{j} = cos(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = cos(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = cos(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);


j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-2.1*pi-1.5*pi,1.9*pi/2));
tic; Z{j} = cos(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = cos(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = cos(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

end
%% tan(x)
if(funcflags.all||funcflags.tan)
    
j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(pi/6,4*pi/9));
tic; Z{j} = tan(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = tan(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = tan(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-3*pi/7,-pi/12));
tic; Z{j} = tan(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = tan(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = tan(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-3*pi/7,4*pi/9));
tic; Z{j} = tan(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = tan(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = tan(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

end
%% sec(x)
if(funcflags.all||funcflags.sec)
    
j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(pi/6,4*pi/9));
tic; Z{j} = sec(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = sec(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sec(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-3*pi/7,-pi/12));
tic; Z{j} = sec(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = sec(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sec(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-3*pi/7,4*pi/9));
tic; Z{j} = sec(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.01:sup(X{j}.x);
zsample{j} = sec(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sec(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

end
%% sqrt(x)
if(funcflags.all||funcflags.sqr)
    
j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(100,10000));
tic; Z{j} = sqrt(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = sqrt(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sqrt(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(0.001,100));
tic; Z{j} = sqrt(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = sqrt(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = sqrt(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

end
%% abs(x)
if(funcflags.all||funcflags.abs)
    
j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(2,4));
tic; Z{j} = abs(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = abs(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);



j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-4,-1));
tic; Z{j} = abs(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = abs(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);


j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(-4,5));
tic; Z{j} = abs(X{j}); tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = abs(xsample{j});
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = abs(x)')
xlabel('x')
ylabel('z')
grid(plots.grid);

end
%% x^a, a positive non-integer
if(funcflags.all||funcflags.pni)
    
j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(0.1,4));
tic; Z{j} = X{j}^0.28; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^0.28;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^{0.28}')
xlabel('x')
ylabel('z')
grid(plots.grid);

j=j+1;
Polyrelax.clear;
X{j} = Polyrelax(Interval(0.1,4));
tic; Z{j} = X{j}^3.7; tocarr(j) = toc;
xsample{j} = inf(X{j}.x):0.1:sup(X{j}.x);
zsample{j} = xsample{j}.^3.7;
Hrep{j} = Polyrelax.Hrep;
P{j} = Polyhedron('A',Hrep{j}.H,'b',Hrep{j}.k,'Ae',Hrep{j}.A,'be',Hrep{j}.b);
figure
plot(xsample{j},zsample{j},'LineWidth',2)
hold on
plot(P{j}.projection([X{j}.i,Z{j}.i]),'Color',plots.color,'Alpha',0.1);
title('z = x^{3.7}')
xlabel('x')
ylabel('z')
grid(plots.grid);

end