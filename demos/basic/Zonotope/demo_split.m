demoID = 'Splitting a zonotope into multiple zonotopes';
util.demoinit;

c = 0.5*ones(2,1);
G = [ 2.5000,    1.2500,    1.5000;
     -0.7500,    0.5000,    3.5000];

X0 = Zonotope(c,G);


fig1 = figure;
xlabel('x_1')
ylabel('x_2')


X0split = split(X0,2);
X0split2 = split(X0,[1,2]);
X0split3 = split(X0,1:3);
X0splitS = split(X0,[1,1,1,1:3]);


figure(fig1);
hold on;
plot(X0,'green',1); 
plot(X0split{1},'yellow',0.3);
plot(X0split{2},'yellow',0.3);

figure;
hold on;
plot(X0,'green',1);
plot(X0split2{1},'yellow',0.3);
plot(X0split2{2},'yellow',0.3);
plot(X0split2{3},'yellow',0.3); 
plot(X0split2{4},'yellow',0.3);


figure;
hold on;
plot(X0,'green',1);
for j=1:size(X0split3,1)
    plot(X0split3{j},'yellow',0.3); 
end

figure;
hold on;
plot(X0,'green',1);
for j=1:size(X0splitS,1)
    plot(X0splitS{j},'yellow',0.3); 
end


