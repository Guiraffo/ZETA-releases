demoID = 'State estimation of descriptor systems using CZs and LZs.';
util.demoinit;

rng('default');
rng(51);



%% System matrices

E = [1, 0, 0;
     0, 1, 0;
     0, 0, 0];
 
A = [0.5,    0,   0;
     0.8, 0.95,   0;
     0.3,  0.1, 0.1];
  
Bu = [1, 0;
     0, 1;
     0, 0];
 
Bw = [0.1,   0,   0;
       0, 1.5,   0;
       0,   0, 0.6];

C = [1,  0, 1;
     1, -1, 0];
 
Dv = [0.5,   0;
       0, 1.5];
   
x0 = [0.5; 0.5; 0.25];   

Nsim = 100;

Ts = 10*pi/Nsim;
simtime = 0:Ts:Ts*Nsim;
simsteps = 0:1:Nsim;

u = [0.5*sin(simtime)+1; 
      -2*cos(simtime)];
   
ngmax = 15;
ncmax = 5;  
   
descsys = DTSystem('descriptor','E',E,'A',A,'Bu',Bu,'Bw',Bw,'C',C,'Dv',Dv);


%% Uncertainties

w_rad = 1;
v_rad = 1;

W = CZonotope(zeros(3,1), w_rad*eye(3));
V = CZonotope(zeros(2,1), v_rad*eye(2));

X0 = CZonotope(x0, diag([0.1, 1.5, 0.6]));
X0_R = LZonotope.realset(3); % For line zonotopes, unbounded initial set

if(~isinside(X0,x0))
    error('The initial condition x0 is unfeasible for this system');
end

%% Settings

% OPTIONS
OPTIONS_CZ_start = util.estimopts('descriptor','u',ngmax,ncmax,[]);
OPTIONS_CZ = util.estimopts('descriptor','pu',ngmax,ncmax,[]);
OPTIONS_LZ_start = util.estimopts('descriptor','u',ngmax,ncmax,[]);
OPTIONS_LZ = util.estimopts('descriptor','pu',ngmax,ncmax,[]);

   
%% Simulate the system

descsys = simulate(descsys,x0,Nsim,W,V,u);
plot(descsys);

   
% Get simulation data
x = descsys.simdata.x;
y = descsys.simdata.y;
u = descsys.simdata.u;

%% Initialize variables

X_CZ = cell(Nsim+1,1);
X_LZ = cell(Nsim+1,1);
Xadm = CZonotope(zeros(3,1), 50*eye(3));

% Initial update
X_CZ{1} = czon_descriptor_estimator(descsys,X0,W,V,[],[],y(:,1),Xadm,OPTIONS_CZ_start);
X_LZ{1} = lzon_descriptor_estimator(descsys,X0_R,LZonotope(W),LZonotope(V),[],[],y(:,1),OPTIONS_LZ_start);

%% State estimation loop
wb = waitbar(0,['Estimation progress: k = ',sprintf('\t\t\t'), num2str(0), '/', num2str(Nsim)]);
for k=2:Nsim+1
    [X_CZ{k},~,timesCZ(k-1)] =  czon_descriptor_estimator(descsys,X_CZ{k-1},W,V,u(:,k),u(:,k-1),y(:,k),Xadm,OPTIONS_CZ);
    [X_LZ{k},~,timesLZ(k-1)] =  lzon_descriptor_estimator(descsys,X_LZ{k-1},LZonotope(W),LZonotope(V),u(:,k),u(:,k-1),y(:,k),OPTIONS_LZ);
    
    % Prints estimation progress each 5 time steps
    if(rem(k-1,5)==0)
        waitbar(k/(Nsim+1),wb,['Estimation progress: k = ',sprintf('\t\t\t'), num2str(k-1), '/', num2str(Nsim)]);
    end        
end
waitbar(1,wb,'Generating figures, please wait.');

%% Radius and interval bounds

RadiusCZ = zeros(1,Nsim+1);
RadiusCZ_1 = zeros(1,Nsim+1);
X1box_CZ = cell(Nsim+1,1);
RadiusLZ = zeros(1,Nsim+1);
RadiusLZ_1 = zeros(1,Nsim+1);
X1box_LZ = cell(Nsim+1,1);

for k=1:Nsim+1
   [RadiusCZ(1,k),RadiusCZ_1(1,k),X1box_CZ{k}] = radius(X_CZ{k});
   [RadiusLZ(1,k),RadiusLZ_1(1,k),X1box_LZ{k}] = radius(X_LZ{k});
end

lowerbounds_CZ = inf([X1box_CZ{:}]);
upperbounds_CZ = sup([X1box_CZ{:}]);
lowerbounds_LZ = inf([X1box_LZ{:}]);
upperbounds_LZ = sup([X1box_LZ{:}]);

%% Plot radius and interval bounds

fig9 = figure;
hold on;
plot(simsteps,RadiusCZ,'b+');
plot(simsteps,RadiusLZ,'rx');
legend('CZ','LZ');
xlabel('k');
ylabel('Radius');
box on;

fig10 = figure;
hold on;
plot(simsteps,RadiusCZ_1,'b+');
plot(simsteps,RadiusLZ_1,'rx');
legend('CZ','LZ');
xlabel('k');
ylabel('Radius: 1-norm');
box on;

figure
hold on;
plot(simsteps,x(1,:),'k');
plot(simsteps,lowerbounds_CZ(1,:),'b--');
plot(simsteps,lowerbounds_LZ(1,:),'r--');
plot(simsteps,upperbounds_CZ(1,:),'b--');
plot(simsteps,upperbounds_LZ(1,:),'r--');
legend('x_1','CZ','LZ');
xlabel('k');
ylabel('x_1');
box on;

figure
hold on;
plot(simsteps,x(2,:),'k');
plot(simsteps,lowerbounds_CZ(2,:),'b--');
plot(simsteps,lowerbounds_LZ(2,:),'r--');
plot(simsteps,upperbounds_CZ(2,:),'b--');
plot(simsteps,upperbounds_LZ(2,:),'r--');
legend('x_2','CZ','LZ');
xlabel('k');
ylabel('x_2');
box on;

figure
hold on;
plot(simsteps,x(3,:),'k');
plot(simsteps,lowerbounds_CZ(3,:),'b--');
plot(simsteps,lowerbounds_LZ(3,:),'r--');
plot(simsteps,upperbounds_CZ(3,:),'b--');
plot(simsteps,upperbounds_LZ(3,:),'r--');
legend('x_3','CZ','LZ');
xlabel('k');
ylabel('x_3');
box on;


%% Computational times

AVGtimes_CZ = mean([timesCZ.total])
AVGtimes_LZ = mean([timesLZ.total])

pause(0.1); close(wb);