function [sys,x,y,w,v] = simulate(sys,x0,varargin)
%SIMULATE simulates a discrete-time system with bounded uncertainties
%
%   SYNTAX: [x,y,w,v] = SIMULATE(sys,x0)
%           [x,y,w,v] = SIMULATE(sys,x0,N)
%           [x,y,w,v] = SIMULATE(sys,x0,N,W,V)
%           [x,y,w,v] = SIMULATE(sys,x0,N,W,V,u)
%
%   INPUTS
%         sys: discrete-time system object (dtsystem)
%          x0: initial state
%           W: bounds on process uncertainties (interval, zonotope, or 
%              czonotope)(optional)
%           V: bounds on measurement uncertainties (interval, zonotope, or 
%              czonotope)(optional)
%           u: known system input (optional)
%           N: the number of time steps to simulate recursively (optional,
%              default = 1 step)
%
%   OUTPUT
%           x: system states
%           y: system measurement
%           w: process uncertainties
%           v: measurement uncertainties

% (C) Copyright 2025 ZETA Developers
%
% This file is a part of the ZETA toolbox
%
%     ZETA is free software: you can redistribute it and/or modify
%     it under the terms of the GNU Lesser General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     ZETA is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU Lesser General Public License for more details.
%
%     You should have received a copy of the GNU Lesser General Public License
%     along with ZETA.  If not, see <http://www.gnu.org/licenses/>.
%
% ZETA on github: https://github.com/Guiraffo/ZETA-releases
%
% Corresponding author: Brenner Santana Rego, brennersr7@usp.br

switch nargin
    case {0,1}
        error('Not enough input arguments.');
    case 2
        N = 1;        
        W = [];
        V = [];        
        u = zeros(0,N+1);
    case 3
        N = varargin{1};
        W = [];
        V = [];        
        u = zeros(0,N+1);        
    case 5
        N = varargin{1};
        W = varargin{2};
        V = varargin{3};
        u = zeros(0,N+1);
    case 6
        N = varargin{1};
        W = varargin{2};
        V = varargin{3};        
        u = varargin{4};
    otherwise
        error_invalid;
end
    
% Struct to reduce MATLAB overhead

dtsys = struct(sys);

% Check dimensions
% Check N
if(~isnumeric(N))
    error_invalid;
elseif(N<1)
    error_invalid;
elseif((size(x0,1)~=dtsys.nx)||(size(x0,2)~=1)) % check x0
    error_dimensions;
elseif((size(u,1)~=dtsys.nu)||(size(u,2)~=(N+1))) % check u and N
    error_dimensions;
else
    % Check W
    isemptyW = isempty(W);
    if(~isemptyW)
        if(isa(W,'Interval'))
            if((size(W,1)~=dtsys.nw)||(size(W,2)~=1))
                error_dimensions;
            end
        elseif(W.dim~=dtsys.nw)
            error_dimensions;
        end
    else
        if(dtsys.nw~=0)
            error_invalid;
        end
    end
    % Check V
    isemptyV = isempty(V);
    if(~isemptyV)
        if(isa(V,'Interval'))
            if((size(V,1)~=dtsys.nv)||(size(V,2)~=1))
                error_dimensions;
            end
        elseif(V.dim~=dtsys.nv)
            error_dimensions;
        end
    else
        if(dtsys.nv~=0)
            error_invalid;
        end
    end    
end

% Generate disturbance samples, if any
if(isemptyW)
    w = zeros(0,N+1);
else
    w = sample(W,N+1);
end

if(isemptyV)
    v = zeros(0,N+1);
else
    v = sample(V,N+1);
end

switch dtsys.subclass
    case 'linear'
        
        % Initialize x and y
        x = zeros(dtsys.nx,N+1); x(:,1) = x0;
        y = zeros(dtsys.ny,N+1); y(:,1) = dtsys.C*x(:,1) + dtsys.Dv*v(:,1);

        % For loop        
        for k=2:N+1
            x(:,k) = dtsys.A*x(:,k-1) + dtsys.Bw*w(:,k-1) + dtsys.Bu*u(:,k-1);
            y(:,k) = dtsys.C*x(:,k)   + dtsys.Dv*v(:,k);
        end
        
    case 'nonlinear'
        
        % Initialize x and y
        x = zeros(dtsys.nx,N+1); x(:,1) = x0;
        y = zeros(dtsys.ny,N+1); y(:,1) = dtsys.g.eval(x(:,1),v(:,1));

        % For loop        
        for k=2:N+1
            x(:,k) = dtsys.f.eval(x(:,k-1),w(:,k-1),u(:,k-1),dtsys.Ts);
            y(:,k) = dtsys.g.eval(x(:,k),v(:,k));
        end        
        
    case 'descriptor'
        
        % Get feasible state and disturbance trajectories
        [x,w] = simulate_descriptor(x0,dtsys.E,dtsys.A,dtsys.Bu,dtsys.Bw,u,w,W,N);
        
        % Initialize y
        y = zeros(dtsys.ny,N+1);        

        % Generate measurements       
        for k=1:N+1
            y(:,k) = dtsys.C*x(:,k) + dtsys.Dv*v(:,k);
        end                
        
    otherwise
        error('DTSystem/simulate is only implemented for the linear and nonlinear cases for now.')
end

sys.simdata = struct('x',x,'y',y,'w',w,'v',v,'u',u);

end

%% Generate feasible trajectories and disturbances of a descriptor system
function [x_seq,w_seq] = simulate_descriptor(x0,E,A,B,Bw,u_seq,w_seed,W,N)
%SIMULATE_DESCRIPTOR generates feasible trajectories and disturbances of a
%                    descriptor system. Requires Yalmip.

W = CZonotope(W);

% Yalmip options
persistent OPTIONS; % Helps computational time
if(isempty(OPTIONS))
    OPTIONS = sdpsettings; % Somehow too heavy to declare in every call
    OPTIONS.verbose = 0;
    OPTIONS.solver = toolsettings.solverlp;
end


nof_x = size(x0,1);
nof_w = size(W.c,1);
nof_gw = size(W.G,2);

x = sdpvar(nof_x,N+1,'full');
csi_w = sdpvar(nof_gw,N+1,'full');


% Constraints on csi_w
%Constraints = norm(reshape(csi_w,nof_gw*(N+1),1),'inf') <= 1;
Constraints = abs(csi_w) <= 1;
%Constraints = abs(csi_w) < 1 + 1e-3;
%Constraints = max(max(abs(csi_w)))<=1;
if(~isempty(W.A)) % If there are equality constraints
    for k=1:N+1
        Constraints = [Constraints;
                       W.A*csi_w(:,k) == W.b];
    end
end

% Constraints on x
Constraints = [Constraints;
               E*x(:,1) == A*x0 + B*u_seq(:,1) + Bw*(W.c+W.G*csi_w(:,1))];               
for k=2:N+1
    Constraints = [Constraints;
                   E*x(:,k) == A*x(:,k-1) + B*u_seq(:,k) + Bw*(W.c+W.G*csi_w(:,k))];
end               

diagnostics = optimize(Constraints,norm(reshape(repmat(W.c,1,N+1)+W.G*csi_w-w_seed,nof_w*(N+1),1),'inf'),OPTIONS);

if(diagnostics.problem~=0) % Yalmip
    error(strcat(['Error in simulate: ',diagnostics.info]));
end  

x_seq = [x0, value(x(:,1:end-1))];
w_seq = repmat(W.c,1,N+1) + W.G*value(csi_w(:,1:end));
               

end


%% Error messaages
function error_invalid
    error('Invalid inputs in dtsystem/simulate.');
end
function error_dimensions
    error('Wrong input dimensions in dtsystem/simulate.');
end

