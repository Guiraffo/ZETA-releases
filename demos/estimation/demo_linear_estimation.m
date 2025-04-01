demoID = 'Linear state estimation using zonotopes and CZs';
util.demoinit;

%% Inputs

rng(30);

X0 = Zonotope([-1;1],[2,   4, 2;
                      2,   0, 2]);
                   
W = Zonotope(Interval(-0.1*ones(2,1),0.1*ones(2,1)));
V = Zonotope(Interval(-0.3*ones(2,1),0.3*ones(2,1)));

linearsys = DTSystem('linear','A',[0, -1;
                                   1,  1],...
                     'Bw',eye(2),...
                     'C',[1, -1;
                          0,  1],'Dv',eye(2));
                               
x0 = [1;-1];

Nsim = 30;

linearsys = simulate(linearsys,x0,Nsim,W,V);

% Initialize cell arrays
X_zon = cell(Nsim+1,1);
X_czon = cell(Nsim+1,1);


% Get simulation data
x = linearsys.simdata.x;
y = linearsys.simdata.y;
u = linearsys.simdata.u;

% Maximum generators and constraints
ngmax = 20;
ncmax = 5;

% OPTIONS
OPTIONS_zon_start = util.estimopts('linear','u',ngmax,ncmax,[]);
OPTIONS_zon = util.estimopts('linear','pu',ngmax,ncmax,[]);
OPTIONS_czon_start = util.estimopts('linear','u',ngmax,ncmax,[]);
OPTIONS_czon = util.estimopts('linear','pu',ngmax,ncmax,[]);

% Initial update step
X_zon{1}  =  zon_linear_estimator(linearsys,X0,W,V,[],y(:,1),OPTIONS_zon_start);
X_czon{1} = czon_linear_estimator(linearsys,CZonotope(X0),CZonotope(W),CZonotope(V),[],y(:,1),OPTIONS_czon_start);


%% For loop
wb = waitbar(0,['Estimation progress: k = ',sprintf('\t\t\t'), num2str(0), '/', num2str(Nsim)]);
for k=2:Nsim+1
    [X_zon{k}, ~,timesZon(k-1)]  =  zon_linear_estimator(linearsys,X_zon{k-1},W,V,u(:,k-1),y(:,k),OPTIONS_zon);
    [X_czon{k},~,timesCZon(k-1)] = czon_linear_estimator(linearsys,X_czon{k-1},CZonotope(W),CZonotope(V),u(:,k-1),y(:,k),OPTIONS_czon);
    
    % Prints estimation progress each 5 time steps
    if(rem(k-1,5)==0)
        waitbar(k/(Nsim+1),wb,['Estimation progress: k = ',sprintf('\t\t\t'), num2str(k-1), '/', num2str(Nsim)]);
    end      
end
waitbar(1,wb,'Generating figures, please wait.');
plot(linearsys);

%% Figures
figure;
hold on;
for k=1:6
    plot(X_zon{k},'red',0.8);
    plot(X_czon{k},'blue',0.8);
end
plot(x(1,1:6),x(2,1:6),'k--x');
legend('Zon','CZ')
xlabel('x_1');
ylabel('x_2');
grid off;
box on;

for k=1:Nsim+1
    [X_zon_rad(k),X_zon_rad1(k),X_zon_hull(:,k)] = radius(X_zon{k});
    [X_czon_rad(k),X_czon_rad1(k),X_czon_hull(:,k)] = radius(X_czon{k});
end

figure;
hold on;
plot(0:1:Nsim,X_zon_rad,'r');
plot(0:1:Nsim,X_czon_rad,'b');
xlabel('k');
ylabel('radius');
legend('Zon','CZ');
box on;

figure;
hold on;
plot(0:1:Nsim,X_zon_rad1,'r');
plot(0:1:Nsim,X_czon_rad1,'b');
xlabel('k');
ylabel('radius: 1-norm');
legend('Zon','CZ');
box on;

LB_X_zon_hull = inf(X_zon_hull);
UB_X_zon_hull = sup(X_zon_hull);
LB_X_czon_hull = inf(X_czon_hull);
UB_X_czon_hull = sup(X_czon_hull);

for j=1:linearsys.nx
    figure;
    hold on;
    plot(0:1:Nsim,LB_X_zon_hull(j,:),'r');
    plot(0:1:Nsim,LB_X_czon_hull(j,:),'b');
    plot(0:1:Nsim,UB_X_zon_hull(j,:),'r');
    plot(0:1:Nsim,UB_X_czon_hull(j,:),'b'); 
    plot(0:1:Nsim,x(j,:),'k--x');
    xlabel('k');
    ylabel(['Interval hull: x_',num2str(j)]);
    legend('Zon','CZ');   
    box on;
end

pause(0.1); close(wb);
