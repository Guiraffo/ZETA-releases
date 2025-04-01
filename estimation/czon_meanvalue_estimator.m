function [Xhat,Xbar,cputimes] = czon_meanvalue_estimator(dtsys,Xhat_prev,W,V,u,y,OPTIONS)
%CZON_MEANVALUE_ESTIMATOR implements a state estimator for nonlinear
%                         discrete-time systems using constrained
%                         zonotopes based on the mean value extension                      

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

if(~(flags.pred || flags.updt))
    error('Invalid estimator operation mode.');
end

% Prediction step
tic;
if(flags.pred)
    Xbar = prediction(dtsys.f,Xhat_prev,W,u,dtsys.Ts);
else
    Xbar = Xhat_prev;
end
cputimes.pred = toc;

% Update step
tic;
if(flags.updt)
    Xhat = update(dtsys.g,Xbar,V,y);
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

function Xbar = prediction(f,X,W,u,Ts)
% Prediction step based on mean value extension using constrained zonotopes (Rego et al., 2021)

% Interval hull
Xhull = intervalhull(X);
if(isempty(W))
    W = zonotope;
    Whull = [];
else
    Whull = intervalhull(W);
end
nof_states = X.dim;

% Jacobian
f_Jacob = f.Jacob(Xhull,Whull,u,Ts);
f_Jacob_x = Interval(f_Jacob(:,1:nof_states)); % Interval call just in case the result is not an interval, so mid and diam still work

% Get approximation point hx
[hx,~,Xbar] = chooseh(X,Xhull,f_Jacob_x,'mindiamm');


% ---------- CZ Zq \supset f(cx,W)  -----------

if(isempty(W.c))
    Zw = f.eval(hx,W.c,u,Ts);
else
    
    f_Jacob_w = f.Jacob(hx,Whull,u,Ts);
    f_Jacob_w = Interval(f_Jacob_w(:,nof_states+1:end)); % Interval call just in case the result is not an interval, so mid and diam still work
    
    % Get approximation point hw
    [hw,~,Wbar] = chooseh(W,Whull,f_Jacob_w,'mindiamm');
    
    % CZ inclusion
    Zwinc = inclusion(CZonotope(W.c-hw,W.G,W.A,W.b),f_Jacob_w,Wbar);
    
    % CZ Zq
    Zw = f.eval(hx,hw,u,Ts) + Zwinc;
    
end

% -----------------------------------------------------------------------

% ----------------------------- CZ Xbar --------------------------------

Zinc = inclusion(CZonotope(X.c-hx,X.G,X.A,X.b),f_Jacob_x,Xbar);

Xbar = Zw + Zinc;

% -----------------------------------------------------------------------


end

function Xhat = update(g,X,V,y)
% Update step based on mean value extension using constrained zonotopes (Rego et al., 2021)

% Interval hull
Xhull = intervalhull(X);
if(isempty(V))
    V = zonotope;
    Vhull = [];
else
    Vhull = intervalhull(V);
end
nof_states = X.dim;

% Jacobian
g_Jacob = g.Jacob(Xhull,Vhull);
g_Jacob_x = Interval(g_Jacob(:,1:nof_states)); % Interval call just in case the result is not an interval, so mid and diam still work

% Get Jtilde
Jtilde = mid(g_Jacob_x);

% Choose h 
[hx,~,Xbar] = chooseh(X,Xhull,Jtilde-g_Jacob_x,'mindiamm');

% ---------- Zonotope Zv \supset lambda(h,V)  -----------

% Affine in V

if(isempty(V))
    Zv = -g.eval(hx,V.c,Ts);
else
    
    % Mean value extension for V
      
    % Jacobian
    g_Jacob_v = g.Jacob(hx,Vhull);
    g_Jacob_v = Interval(g_Jacob_v(:,nof_states+1:end)); % Interval call just in case the result is not an interval, so mid and diam still work
           
    % Choose h 
    [hv,~,Vbar] = chooseh(V,Vhull,g_Jacob_v,'mindiamm');    
    Zvinc = inclusion(CZonotope(V.c-hv,V.G,V.A,V.b),g_Jacob_v,Vbar);
    
    % It actually must contain minus the set
    Zv = -(g.eval(hx,hv) + Zvinc);
    
    
end

% -----------------------------------------------------------------------

% ---------------------- Constrained zonotope Y -------------------------

Zinc = inclusion(CZonotope(X.c-hx,X.G,X.A,X.b),(Jtilde-g_Jacob_x),Xbar);

Y = y + Jtilde*hx + Zv + Zinc;

% -----------------------------------------------------------------------

% Generalized intersection

Xhat = intersection(X,Y,Jtilde);

end



