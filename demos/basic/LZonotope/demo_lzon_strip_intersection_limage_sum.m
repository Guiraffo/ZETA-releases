demoID = 'Intersection, affine transformation, and Minkowsky sum of line zonotopes';
util.demoinit;

%% Inputs

S1 = Strip([1; -1; 1], 1, 0.1);
S2 = Strip([1; 1; 1], 1, 0.1);

%% Intersection

S1_LZ = LZonotope(S1);
S2_LZ = LZonotope(S2);
S1S2inter = intersection(S1_LZ,S2_LZ);

%% Linear transformation

R = [1, 0, 1; 0, 1, 1; 0, 1, -1];
q = -ones(3,1);
S3 = R*S1S2inter + q;

%% Minkowski sum

S4 = S1S2inter + S3;

%% Figures 

figure;
hold on;
axis0 = [-0.1900, 2.1900, -0.7281, 0.7281, -0.3960, 0.3960];
axis(axis0);
view([-37.5,30]);
plot(S1,'blue',0.3);
plot(S2,'red',0.3);
axis(axis0);
plot(S1S2inter,'black',1);
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
legend('S_1','S_2','S_1 \cap S_2 (CLG-rep)')

figure;
hold on;
axis0 = 1*[-1.1900, 2.1900, -2.7281, 1.7281, -1.3960, 1.3960];
axis(axis0);
view([-37.5,30]);
plot(S1S2inter,'black',1);
plot(S3,'yellow',1);
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
legend('S_1 \cap S_2','S_3')
title('S_3 = R(S_1 \cap S_2) + q');

figure;
hold on;
axis0 = 1*[-1.1900, 2.1900, -4.7281, 1.7281, -1.3960, 1.3960];
axis(axis0);
view([-37.5,30]);
plot(S4,'red',0.8);
xlabel('x_1');
ylabel('x_2');
zlabel('x_3');
legend('S_4')
title('S_4 = (S_1 \cap S_2) + S_3')
