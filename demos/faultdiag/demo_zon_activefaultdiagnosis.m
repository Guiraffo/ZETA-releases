demoID = 'Active fault diagnosis of linear systems using zonotopes';
util.demoinit;

%% Inputs

nof_inputs = 2;
u_rad = 9;

[U.H,U.k] = hrep(Interval(-u_rad*ones(nof_inputs,1),u_rad*ones(nof_inputs,1)));

A = cell(5,1);
Bu = cell(5,1);
Bw = cell(5,1);
C = cell(5,1);
Dv = cell(5,1);
r = cell(5,1);
s = cell(5,1);

A{1} = [0.6, 0.2; -0.4, -0.2];
A{2} = A{1};
A{3} = A{1};
A{4} = [1.2, 0.2; -0.4, -0.2];
A{5} = [2.0, 0.2; -0.4, -0.7];

Bu{1} = eye(2);
Bu{2} = [1, 0; 0, 0];
Bu{3} = [0, 0; 0, 1];
Bu{4} = Bu{1};
Bu{5} = Bu{1};

Bw{1} = [1, 0; 0, 0];
Bw{2} = Bw{1};
Bw{3} = Bw{1};
Bw{4} = Bw{1};
Bw{5} = Bw{1};

C{1} = [1, 0];
C{2} = C{1};
C{3} = C{1};
C{4} = C{1};
C{5} = C{1};

Dv{1} = 1;
Dv{2} = Dv{1};
Dv{3} = Dv{1};
Dv{4} = Dv{1};
Dv{5} = Dv{1};

r{1} = zeros(2,1);
r{2} = r{1};
r{3} = r{1};
r{4} = r{1};
r{5} = r{1};

s{1} = 0;
s{2} = s{1};
s{3} = s{1};
s{4} = s{1};
s{5} = s{1};

for j=1:5
    ltisystems(j) = DTSystem('linear','A',A{j},'Bu',Bu{j},'Bw',Bw{j},'C',C{j},'Dv',Dv{j});
end


X0 = Zonotope(-3*ones(2,1), 0.2*eye(2));
W = Zonotope(zeros(2,1), 0.5*eye(2));
V = Zonotope(0, 0.2);

N = 2;

%% Design of separating inputs

OPTIONS_AFD.quad = util.faultdiagopts('quadratic',eye(nof_inputs),N,2,1.0e-4);
OPTIONS_AFD.lin = util.faultdiagopts('linear',eye(nof_inputs),N,2,1.0e-4);

[usep_quad, deltahat_m_quad, cputimes_quad] = zonAFD_inputdesign(ltisystems,r,s,X0,W,V,U,OPTIONS_AFD.quad);
[usep_lin,  deltahat_m_lin , cputimes_lin]  = zonAFD_inputdesign(ltisystems,r,s,X0,W,V,U,OPTIONS_AFD.lin);

%% Compute reachable sets

Psi_quad = cell(5,1);
Psi_lin = cell(5,1);
useparating_quad = reshape(usep_quad,nof_inputs,N);
useparating_lin =  reshape(usep_lin,nof_inputs,N);

for j=1:5
    [~,Psi_quad{j},~] = AFD.LTI_reachable_sets(r{j},A{j},Bu{j},Bw{j},s{j},C{j},Dv{j},X0,W,V,useparating_quad);
    [~,Psi_lin{j},~]  = AFD.LTI_reachable_sets(r{j},A{j},Bu{j},Bw{j},s{j},C{j},Dv{j},X0,W,V,useparating_lin);
end

%% Verify if intersection is empty

nof_models = 5;

Zinter_quad = cell(nof_models,nof_models);
Zinter_lin = cell(nof_models,nof_models);
flagintersection_quad = cell(nof_models,nof_models);
flagintersection_lin = cell(nof_models,nof_models);

for i=1:nof_models-1
    for j=i+1:nof_models
        Zinter_quad{i,j} = intersection(CZonotope(Psi_quad{i}),CZonotope(Psi_quad{j}));
        flagintersection_quad{i,j} = isempty(Zinter_quad{i,j});
        Zinter_lin{i,j} = intersection(CZonotope(Psi_lin{i}),CZonotope(Psi_lin{j}));
        flagintersection_lin{i,j} = isempty(Zinter_lin{i,j});        
    end
end

flagintersection_quad
cputimes_quad
deltahat_m_quad
separatinginput_quad = usep_quad
norminput_quad = norm(usep_quad,2)

flagintersection_lin
cputimes_lin
deltahat_m_lin
separatinginput_lin = usep_lin
norminput_lin = norm(usep_lin,1)

%% Figures

figure;
hold on;
colorarr = {'green','cyan','yellow','red','black'};
legendarr = cell(length(ltisystems),1);
for j=1:length(ltisystems)
    plot(Psi_quad{j},colorarr{j},1);
    legendarr{j} = ['Model ',num2str(j)];
end
legend(legendarr{:},'Location','southeast');
box on;
grid off;
xlabel('y_1');
ylabel('y_2');
title('Separated output reachable tubes, MIQP');
axis([-10 -1 -25 2])

figure;
hold on;
colorarr = {'green','cyan','yellow','red','black'};
legendarr = cell(length(ltisystems),1);
for j=1:length(ltisystems)
    plot(Psi_lin{j},colorarr{j},1);
    legendarr{j} = ['Model ',num2str(j)];
end
legend(legendarr{:},'Location','southeast');
box on;
grid off;
xlabel('y_1');
ylabel('y_2');
title('Separated output reachable tubes, MILP');
axis([-10 -1 -25 2])    