demoID = 'Linear image, Minkowski sum, and generalized intersection of constrained zonotopes';
util.demoinit;


%% Inputs

Z = CZonotope([0.1763; 0.7954; 0.7831], [ 0.6317, -0.2426, -0.6123,  0.5660, -0.9267;
                                         -0.9282,  0.0370, -0.4554,  0.7007, -0.7666;
                                          0.3835,  0.3159,  0.4372,  0.5505,  0.5026], [-0.5216, -0.4904, 0.7153, 0.8996, 0.1234], -0.6424);
                                      
W = CZonotope([0.5405;-0.0152; 0.2625], [ 0.6790,  0.3588, -0.8654, -0.3416;
                                         -0.0779,  0.3016,  0.5429,  0.0213;
                                         -0.0041, -0.4624, -0.0380, -0.4727], [-0.3790, 0.2537, 0.1149, -0.3628], -0.2103);
                                      
Y = CZonotope([ -0.4841; 0.1645], [ -0.6767,  0.6516,  0.4686,  0.5574;
                                     0.1963, -0.6872, -0.1827,  0.6079], [0.5721, 0.1846, 0.3290, 0.2931], -0.1487);
                                     
R = [0.5136, 0.0371, 0.6204;
     0.5013, 0.7081, 0.7778];

%% Operations
 
Zlimage = R*Z;
Zsum = Z+W;
ZintW = Z&W;
ZgenintY = intersection(Z,Y,R);


%% Figures

figure
hold on;
plot(Z,'black',0.3);
plot(W,'green',0.3);
plot(ZintW,'blue',1,'fast');
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')
legend('Z','W','Z \cap W');
title('Intersection of two CZs Z \cap W')

figure
plot(Zlimage,'magenta',0.3);
title('Linear image R*Z')
xlabel('x_1')
ylabel('x_2')


figure
plot(Zsum,'blue',0.3);
title('Minkowski sum Z+W')
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')

figure
plot(Y,'cyan',0.3);
xlabel('x_1')
ylabel('x_2')
legend('Y');
title('Set Y for generalized intersection')

figure
hold on
plot(Z,'black',0.3);
plot(ZgenintY,'yellow',1,'Yalmip');
xlabel('x_1')
ylabel('x_2')
zlabel('x_3')
legend('Z','Z \cap_R Y');
title('Generalized intersection Z \cap_R Y')