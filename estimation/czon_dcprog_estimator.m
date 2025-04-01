function [Xhat,Xbar,cputimes] = czon_dcprog_estimator(dtsys,Xhat_prev,W,V,u,y,OPTIONS)
% Nonlinear estimator using constrained zonotopes based on DC programming as in de Paula et al (2024)

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

% Setup algorithm flags based on input string
flags.pred = any(OPTIONS.alg=='p');
flags.updt = any(OPTIONS.alg=='u');

% Prediction step
tic;
if(flags.pred)
    Xbar = prediction(dtsys.f,Xhat_prev,W,u,dtsys.Ts,OPTIONS);
else
    Xbar = Xhat_prev;
end
cputimes.pred = toc;

% Update step
tic;
if(flags.updt)
    Xhat = update(dtsys.g,Xbar,V,y,OPTIONS);
else
    Xhat = Xbar;
end
cputimes.updt = toc;

% Complexity reduction
tic;
Xhat = reduction(Xhat,OPTIONS.ngmax,OPTIONS.ncmax);
cputimes.redu = toc;

% Total computational time
cputimes.total = cputimes.pred + cputimes.updt + cputimes.redu;

end

%% Prediction step
function [Xbar,hx,hw,P_vert] = prediction(modeldc,X,W,u,Ts,OPTIONS)
% Implementation of the nonlinear prediction step based on DC programming (de Paula et al, 2024)

nof_x = X.dim;
if(isempty(W))
    W = CZonotope;
end
Z = [X;W];
  
switch OPTIONS.vertices
    case 'hull'       
        P = Zonotope(intervalhull(Z));
        P_vert = util.partope_vertices(P);
        hz = P.c;            
    case 'partope'        
        P = partopebound(Z);
        P_vert = util.partope_vertices(P);
        hz = P.c;       
    otherwise
        error('Invalid algorithm option in DC programming propagate.');
end


switch OPTIONS.decomposition
    case 'exact'
        lambda = [];
    case 'aBB'
        lambda = dcprog.getlambda_alphaBB(modeldc,intervalhull(Z),nof_x,u,Ts);
    otherwise
        error('Invalid decomposition mode.');
end
hx = hz(1:nof_x);
hw = hz(nof_x+1:end);

% Get remainder using DC programming
[e_minus,e_plus] = get_error_bounds(modeldc,P_vert,hx,hw,OPTIONS,lambda,u,Ts);
Rem = Zonotope(0.5*(e_plus+e_minus), 0.5*diag(e_plus-e_minus));

f_ath = modeldc.eval(hx,hw,u,Ts);
f_Jacob_ath = modeldc.Jacob(hx,hw,u,Ts);

Xbar = (f_ath - f_Jacob_ath*hz) + f_Jacob_ath*Z + Rem;

end

%% Update step
function [Xhat,hx,hv,P_vert] = update(modeldc,X,V,y,OPTIONS)
% Implementation of the nonlinear update step based on DC programming (de Paula et al, 2024)

nof_x = X.dim;
if(isempty(V))
    V = CZonotope;
end
Z = [X;V];

switch OPTIONS.vertices
    case 'hull'       
        P = Zonotope(intervalhull(Z));
        P_vert = util.partope_vertices(P);
        hz = P.c;            
    case 'partope'        
        P = partopebound(Z);
        P_vert = util.partope_vertices(P);
        hz = P.c;       
    otherwise
        error('Invalid algorithm option in DC programming propagate.');
end
hx = hz(1:nof_x);
hv = hz(nof_x+1:end);

switch OPTIONS.decomposition
    case 'exact'
        lambda = [];
    case 'aBB'
        lambda = dcprog.getlambda_alphaBB(modeldc,intervalhull(Z),nof_x);
    otherwise
        error('Invalid decomposition mode.');
end

% Get remainder using DC programming
[e_minus,e_plus] = get_error_bounds(modeldc,P_vert,hx,hv,OPTIONS,lambda);
Rem = Zonotope(0.5*(e_plus+e_minus), 0.5*diag(e_plus-e_minus));

g_ath = modeldc.eval(hx,hv);
g_Jacob_ath = modeldc.Jacob(hx,hv);
g_Jacob_ath_x = g_Jacob_ath(:,1:nof_x);
g_Jacob_ath_v = g_Jacob_ath(:,nof_x+1:end);

Y = (y - g_ath + g_Jacob_ath*hz) - g_Jacob_ath_v*V - Rem;
Xhat = intersection(X,Y,g_Jacob_ath_x);

end

%% Generate linearization error bounds using DC programming
function [e_minus,e_plus] = get_error_bounds(modeldc,vert,hx,hw,OPTIONS,lambda,varargin)
% Compute bounds on the remainder as in Alamo et al (2008) and de Paula et al (2024)

    nof_vertices = size(vert,2);
    nof_x = size(hx,1);

    err_minor = [];
    err_major = [];
    for j=1:nof_vertices
        x = vert(1:nof_x,j);
        w = vert(nof_x+1:end,j);
        err_minor(:,j) = error_minorant(modeldc,x,w,hx,hw,OPTIONS,lambda,varargin{:});
        err_major(:,j) = error_majorant(modeldc,x,w,hx,hw,OPTIONS,lambda,varargin{:});
    end
    e_minus = min(err_minor,[],2);
    e_plus = max(err_major,[],2);

end

%% Calculate error minorant and majorant
function out = error_minorant(modeldc,x,w,hx,hw,OPTIONS,lambda,varargin)
% Compute the candidate minorant at a vertex as in Alamo et al (2008) and de Paula et al (2024)

linearized_func_at_xw = linearized_function(modeldc,x,w,hx,hw,varargin{:});
switch OPTIONS.decomposition
    case 'exact'
        convxb_at_xw = modeldc.dcB(x,w,varargin{:});
        linearized_convxa_at_xw = linearized_function_convxa(modeldc,x,w,hx,hw,varargin{:});
    case 'aBB'
        convxb_at_xw = evalfunction_convxbABB(modeldc,x,w,lambda,varargin{:});
        linearized_convxa_at_xw = linearized_function_convxaABB(modeldc,x,w,hx,hw,lambda,varargin{:});
    otherwise
        error('Invalid decomposition mode.');
end

out = linearized_convxa_at_xw - convxb_at_xw - linearized_func_at_xw;
    
end
function out = error_majorant(modeldc,x,w,hx,hw,OPTIONS,lambda,varargin)
% Compute the candidate majorant at a vertex as in Alamo et al (2008) and de Paula et al (2024)

linearized_func_at_xw = linearized_function(modeldc,x,w,hx,hw,varargin{:});
switch OPTIONS.decomposition
    case 'exact'
        convxa_at_xw = modeldc.daA(x,w,varargin{:});
        linearized_convxb_at_xw = linearized_function_convxb(modeldc,x,w,hx,hw,varargin{:});
    case 'aBB'
        convxa_at_xw = evalfunction_convxaABB(modeldc,x,w,lambda,varargin{:});
        linearized_convxb_at_xw = linearized_function_convxbABB(modeldc,x,w,hx,hw,lambda,varargin{:});
    otherwise
        error('Invalid decomposition mode.');
end    

out = convxa_at_xw - linearized_convxb_at_xw - linearized_func_at_xw;
    
end

%% Evaluate linearized main function
function out = linearized_function(modeldc,x,w,hx,hw,varargin)
    out = modeldc.eval(hx,hw,varargin{:}) + modeldc.Jacob(hx,hw,varargin{:})*([x;w]-[hx;hw]);
end
%% Evaluation of fa and fb using exact DC decomposition, plus linearized versions
function out = linearized_function_convxa(modeldc,x,w,hx,hw,varargin)
    out = modeldc.dcA(hx,hw,varargin{:}) + modeldc.dcA_Jacob(hx,hw,varargin{:})*([x;w]-[hx;hw]);
end
function out = linearized_function_convxb(modeldc,x,w,hx,hw,varargin)
    out = modeldc.dcB(hx,hw,varargin{:}) + modeldc.dcB_Jacob(hx,hw,varargin{:})*([x;w]-[hx;hw]);
end
%% Evaluation of fa and fb, and respective Jacobians, using the alpha-BB method, plus linearized versions
function out = evalfunction_convxaABB(modeldc,x,w,lambda,varargin)
% fconvxa for system DC decomposition using alphaBB method: fa = f + fb, fb = 0.5*lambda*(z'*z)
z = [x;w];
out = modeldc.eval(x,w,varargin{:}) + 0.5*lambda*(z.'*z);
end
function out = evalfunction_convxbABB(modeldc,x,w,lambda,varargin)
% fconvxb for system DC decomposition using alphaBB convexification: fa = f + fb, fb = 0.5*lambda*(z'*z)
z = [x;w];
out = 0.5*lambda*(z.'*z);
end
function out = evalfunction_convxaABB_Jacob(modeldc,x,w,lambda,varargin)
% fconvxa (Jacobian) for system DC decomposition using alphaBB method: fa = f + fb, fb = 0.5*lambda*(z'*z)
z = [x;w];
out = modeldc.Jacob(x,w,varargin{:}) + lambda*z.';
end
function out = evalfunction_convxbABB_Jacob(modeldc,x,w,lambda,varargin)
% fconvxb (Jacobian) for system DC decomposition using alphaBB method: fa = f + fb, fb = 0.5*lambda*(z'*z)
z = [x;w];
out = lambda*z.';
end
function out = linearized_function_convxaABB(modeldc,x,w,hx,hw,lambda,varargin)
    out = evalfunction_convxaABB(modeldc,hx,hw,lambda,varargin{:}) + evalfunction_convxaABB_Jacob(modeldc,hx,hw,lambda,varargin{:})*([x;w]-[hx;hw]);
end
function out = linearized_function_convxbABB(modeldc,x,w,hx,hw,lambda,varargin)
    out = evalfunction_convxbABB(modeldc,hx,hw,lambda,varargin{:}) + evalfunction_convxbABB_Jacob(modeldc,hx,hw,lambda,varargin{:})*([x;w]-[hx;hw]);
end
