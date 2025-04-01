demoID = 'Nonlinear state estimation using constrained zonotopes';
util.demoinit;

%% Inputs

rng(51);

X0 = CZonotope([5;0.5],[0.5,   1, -0.5;
                        0.5, 0.5,    0]);
                                     
W = CZonotope(zeros(2,1),0.5*eye(2));
V = CZonotope(zeros(2,1),0.2*eye(2));

nonlinearsys = DTSystem('nonlinear','Model1',2,0,2,2,2);
                               
x0 = [5.2;0.65];

Nsim = 100;

nonlinearsys = simulate(nonlinearsys,x0,Nsim,W,V);

%%

% Initialize cell arrays
X_CZMV = cell(Nsim+1,1);
X_CZFO = cell(Nsim+1,1);
X_CZPR = cell(Nsim+1,1);
X_CZDC = cell(Nsim+1,1);

% Get simulation data
x = nonlinearsys.simdata.x;
y = nonlinearsys.simdata.y;
u = nonlinearsys.simdata.u;

% Maximum generators and constraints
ngmax = 20;
ncmax = 8;
ngprior = [];   

% OPTIONS
OPTIONS_CZMV_start = util.estimopts('CZMV','u',ngmax,ncmax,[]);
OPTIONS_CZMV = util.estimopts('CZMV','pu',ngmax,ncmax,[]);
OPTIONS_CZFO_start = util.estimopts('CZFO','u',ngmax,ncmax,ngprior);
OPTIONS_CZFO = util.estimopts('CZFO','pu',ngmax,ncmax,ngprior);
OPTIONS_CZPR_start = util.estimopts('CZPR','u',ngmax,ncmax,[]);
OPTIONS_CZPR = util.estimopts('CZPR','pu',ngmax,ncmax,[]);
OPTIONS_CZDC_start = util.estimopts('CZDC','u',ngmax,ncmax,[]);
OPTIONS_CZDC = util.estimopts('CZDC','pu',ngmax,ncmax,[]);

% Initial update step
X_CZMV{1} =  czon_meanvalue_estimator(nonlinearsys,X0,W,V,[],y(:,1),OPTIONS_CZMV_start);
X_CZFO{1} = czon_firstorder_estimator(nonlinearsys,X0,W,V,[],y(:,1),OPTIONS_CZFO_start);
X_CZPR{1} =  czon_polyrelax_estimator(nonlinearsys,X0,W,V,[],y(:,1),OPTIONS_CZPR_start);
X_CZDC{1} =     czon_dcprog_estimator(nonlinearsys,X0,W,V,[],y(:,1),OPTIONS_CZDC_start);


%% State estimation loop
wb = waitbar(0,['Estimation progress: k = ',sprintf('\t\t\t'), num2str(0), '/', num2str(Nsim)]);
for k=2:Nsim+1
    [X_CZMV{k},~,timesCZMV(k-1)] =  czon_meanvalue_estimator(nonlinearsys,X_CZMV{k-1},W,V,u(:,k-1),y(:,k),OPTIONS_CZMV);
    [X_CZFO{k},~,timesCZFO(k-1)] = czon_firstorder_estimator(nonlinearsys,X_CZFO{k-1},W,V,u(:,k-1),y(:,k),OPTIONS_CZFO);
    [X_CZPR{k},  timesCZPR(k-1)] =  czon_polyrelax_estimator(nonlinearsys,X_CZPR{k-1},W,V,u(:,k-1),y(:,k),OPTIONS_CZPR);
    [X_CZDC{k},~,timesCZDC(k-1)] =     czon_dcprog_estimator(nonlinearsys,X_CZDC{k-1},W,V,u(:,k-1),y(:,k),OPTIONS_CZDC);
    
    % Prints estimation progress each 5 time steps
    if(rem(k-1,5)==0)
        waitbar(k/(Nsim+1),wb,['Estimation progress: k = ',sprintf('\t\t\t'), num2str(k-1), '/', num2str(Nsim)]);
    end    
end
waitbar(1,wb,'Generating figures, please wait.');
plot(nonlinearsys);

%% Figures
figure;
hold on;
for k=1:6
    plot(X_CZFO{k},'red',0.8);
    plot(X_CZDC{k},'cyan',0.8);      
    plot(X_CZMV{k},'blue',0.8);    
    plot(X_CZPR{k},'magenta',0.8);        
end
plot(x(1,1:6),x(2,1:6),'k--x');
legend('CZFO','CZDC','CZMV','CZPR')
xlabel('x_1');
ylabel('x_2');
grid off;
box on;

for k=1:Nsim+1
    [X_CZMV_rad(k),X_CZMV_rad1(k),X_CZMV_hull(:,k)] = radius(X_CZMV{k});
    [X_CZFO_rad(k),X_CZFO_rad1(k),X_CZFO_hull(:,k)] = radius(X_CZFO{k});
    [X_CZPR_rad(k),X_CZPR_rad1(k),X_CZPR_hull(:,k)] = radius(X_CZPR{k});
    [X_CZDC_rad(k),X_CZDC_rad1(k),X_CZDC_hull(:,k)] = radius(X_CZDC{k});
    
    X_CZMV_ivol(k) = volume(X_CZMV{k},'box-nthroot');
    X_CZMV_pvol(k) = volume(X_CZMV{k},'partope-nthroot');
    X_CZFO_ivol(k) = volume(X_CZFO{k},'box-nthroot');
    X_CZFO_pvol(k) = volume(X_CZFO{k},'partope-nthroot');
    X_CZPR_ivol(k) = volume(X_CZPR{k},'box-nthroot');
    X_CZPR_pvol(k) = volume(X_CZPR{k},'partope-nthroot');    
    X_CZDC_ivol(k) = volume(X_CZDC{k},'box-nthroot');
    X_CZDC_pvol(k) = volume(X_CZDC{k},'partope-nthroot');        
end

figure;
hold on;
plot(0:1:Nsim,X_CZMV_rad,'r');
plot(0:1:Nsim,X_CZFO_rad,'b');
plot(0:1:Nsim,X_CZPR_rad,'m');
plot(0:1:Nsim,X_CZDC_rad,'k');
xlabel('k');
ylabel('radius');
legend('CZMV','CZFO','CZPR','CZDC');
box on;

figure;
hold on;
plot(0:1:Nsim,X_CZMV_rad1,'r');
plot(0:1:Nsim,X_CZFO_rad1,'b');
plot(0:1:Nsim,X_CZPR_rad1,'m');
plot(0:1:Nsim,X_CZDC_rad1,'k');
xlabel('k');
ylabel('radius: 1-norm');
legend('CZMV','CZFO','CZPR','CZDC');
box on;

figure;
hold on;
plot(0:1:Nsim,X_CZMV_ivol,'r');
plot(0:1:Nsim,X_CZFO_ivol,'b');
plot(0:1:Nsim,X_CZPR_ivol,'m');
plot(0:1:Nsim,X_CZDC_ivol,'k');
xlabel('k');
ylabel('volume: box-nthroot');
legend('CZMV','CZFO','CZPR','CZDC');
box on;

figure;
hold on;
plot(0:1:Nsim,X_CZMV_pvol,'r');
plot(0:1:Nsim,X_CZFO_pvol,'b');
plot(0:1:Nsim,X_CZPR_pvol,'m');
plot(0:1:Nsim,X_CZDC_pvol,'k');
xlabel('k');
ylabel('volume: partope-nthroot');
legend('CZMV','CZFO','CZPR','CZDC');
box on;


LB_X_CZMV_hull = inf(X_CZMV_hull);
UB_X_CZMV_hull = sup(X_CZMV_hull);
LB_X_CZFO_hull = inf(X_CZFO_hull);
UB_X_CZFO_hull = sup(X_CZFO_hull);
LB_X_CZPR_hull = inf(X_CZPR_hull);
UB_X_CZPR_hull = sup(X_CZPR_hull);
LB_X_CZDC_hull = inf(X_CZDC_hull);
UB_X_CZDC_hull = sup(X_CZDC_hull);

for j=1:nonlinearsys.nx
    figure;
    hold on;
    plot(0:1:Nsim,LB_X_CZMV_hull(j,:),'r');
    plot(0:1:Nsim,LB_X_CZFO_hull(j,:),'b');
    plot(0:1:Nsim,LB_X_CZPR_hull(j,:),'m');
    plot(0:1:Nsim,LB_X_CZDC_hull(j,:),'k');
    plot(0:1:Nsim,UB_X_CZMV_hull(j,:),'r');
    plot(0:1:Nsim,UB_X_CZFO_hull(j,:),'b');
    plot(0:1:Nsim,UB_X_CZPR_hull(j,:),'m');    
    plot(0:1:Nsim,UB_X_CZDC_hull(j,:),'k');    
    plot(0:1:Nsim,x(j,:),'k--x');
    xlabel('k');
    ylabel(['Interval hull: x_',num2str(j)]);
    legend('CZMV','CZFO','CZPR','CZDC');  
    box on;
end

AVGtimes_CZMV = mean([timesCZMV.total])
AVGtimes_CZFO = mean([timesCZFO.total])
AVGtimes_CZPR = mean([timesCZPR.total])
AVGtimes_CZDC = mean([timesCZDC.total])

pause(0.1); close(wb);

