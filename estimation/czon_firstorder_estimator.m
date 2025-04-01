function [Xhat,Xbar,cputimes] = czon_firstorder_estimator(dtsys,Xhat_prev,W,V,u,y,OPTIONS)
%CZON_FIRSTORDER_ESTIMATOR implements a state estimator for nonlinear
%                          discrete-time systems using constrained
%                          zonotopes based on the first-order Taylor
%                          extension

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

% Complexity reduction
tic;
Xbar = reduction(Xbar,OPTIONS.ngmax,OPTIONS.ncmax,OPTIONS.ngprior);
cputimes.redu = toc;

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
Xhat = reduction(Xhat,OPTIONS.ngmax,OPTIONS.ncmax,OPTIONS.ngprior);
cputimes.redu = cputimes.redu + toc;

% Total computational time
cputimes.total = cputimes.pred + cputimes.updt + cputimes.redu;

end

function Xbar = prediction(f,X,W,u,Ts)
% Prediction step based on first-order Taylor extension using constrained zonotopes (Rego et al., 2021)

% Interval hull
Xhull = intervalhull(X);
if(isempty(W))
    W = zonotope;
    Whull = [];
else
    Whull = intervalhull(W);
end
nof_states = X.dim;

Z = [X;W];
Zhull = [Xhull;Whull];

% Choose h
hz = chooseh(Z,Zhull,[],'distcenter');
hx = hz(1:nof_states,1);
hw = hz((nof_states+1):end,1);

% Evaluate the nonlinear function, Jacobian, and triangular Hessians

fh = f.eval(hx,hw,u,Ts);
Jacobfh = f.Jacob(hx,hw,u,Ts);
HessianTriangHull = f.HessianTriang(Xhull,Whull,u,Ts);
    

% Initialize auxiliary variables
in_dim = Z.dim;
in_ng = Z.ng;
in_nc = Z.nc;
out_dim = size(fh,1);
Qbar = cell(out_dim,1);
cR = zeros(out_dim,1);
Gq = zeros(out_dim,sum(1:in_ng));
Gh = zeros(out_dim,out_dim);

% Generators
for q=1:out_dim
   HessianTriangHull{q} = Interval(HessianTriangHull{q}); % Interval call just in case the result is not an interval, so mid and diam still work  
    
   Qbar{q} = Z.G.'*HessianTriangHull{q}*Z.G;
   MidQbar = mid(Qbar{q});
   RadQbar = rad(Qbar{q});
   
   % cR
   cR(q) = (1/2)*trace(MidQbar);
     
   % Gq
   MidQdouble = MidQbar + MidQbar.';
   Diag = (1/4)*diag(MidQdouble).';
  
   UppTriang = zeros(1,sum(1:size(Z.G,2)-1)); counter = 0;
   for j=2:size(Z.G,2)
       for i = 1:(j-1)
           
           counter = counter + 1;
           UppTriang(counter) = MidQdouble(i,j);
       
       end
   end 
   
   Gq(q,:) = [Diag, UppTriang];
   
   % Gh
   Gh(q,q) = sum(sum(abs(RadQbar)));
    
end

% Constraints
Atildezeta = zeros(sum(1:in_nc),in_ng);
Atildecsi = zeros(sum(1:in_nc),sum(1:(in_ng-1)));
btilde = zeros(sum(1:in_nc),1);

counterROW = 0; % Row counter associated to \forall r <= s
for s=1:in_nc
    for r=1:s
        
        counterROW = counterROW + 1;       
        for i=1:in_ng
            Atildezeta(counterROW,i) = (1/2)*Z.A(r,i)*Z.A(s,i);
        end
        
        counterCOL = 0; % Column counter for Atildecsi, associated with \forall i<j
        for j=2:in_ng
            for i=1:(j-1)
                counterCOL = counterCOL + 1;
                Atildecsi(counterROW,counterCOL) = Z.A(r,i)*Z.A(s,j) + Z.A(r,j)*Z.A(s,i);               
            end
        end   
        
        btilde(counterROW) = Z.b(r)*Z.b(s) - sum(Atildezeta(counterROW,:),2); % The last sum take advantage of the summands being the same elements from Atildezeta
            
    end
end

p = Z.c - hz;
Atilde = [Atildezeta, Atildecsi, zeros(sum(1:in_nc), out_dim)];

% Xbar = CZonotope(fh + Jacobfh*p + cR, [Jacobfh*Z.G, Gq, Gh],...
%                   [Z.A, zeros(in_nc,size(Gq,2)+size(Gh,2)); zeros(size(Atilde,1), in_ng), Atilde],...
%                   [Z.b; btilde]); 

Rem = CZonotope(cR, [Gq, Gh], Atilde, btilde);
if(norm(p)>1e-13) % Tolerance value for zero. If p is not zero, add the inclusion part
    L = Interval(zeros(out_dim,in_dim),zeros(out_dim,in_dim));
    for q=1:out_dim
        L(q,:) = p.'*HessianTriangHull{q};
    end    
    Rem = Rem + inclusion(CZonotope(p,2*Z.G,Z.A,Z.b),L);   
end

Xbar = (fh - Jacobfh*hz) + Jacobfh*Z + Rem;

% -----------------------------------------------------------------------

end

function Xhat = update(g,X,V,y)
% Update step based on firt-order Taylor extension using constrained zonotopes (Rego et al., 2021)

% Interval hull
Xhull = intervalhull(X);
if(isempty(V))
    V = zonotope;
    Vhull = [];
else
    Vhull = intervalhull(V);
end
nof_states = X.dim;

Z = [X;V];
Zhull = [Xhull;Vhull];

% Choose h
hz = chooseh(Z,Zhull,[],'distcenter');
hx = hz(1:nof_states,1);
hv = hz((nof_states+1):end,1);

% Evaluate the nonlinear function, Jacobian, and triangular Hessians

gh = g.eval(hx,hv);
Jacobgh = g.Jacob(hx,hv);
HessianTriangHull = g.HessianTriang(Xhull,Vhull);
    

% Initialize auxiliary variables
in_dim = Z.dim;
in_ng = Z.ng;
in_nc = Z.nc;
out_dim = size(gh,1);
Qbar = cell(out_dim,1);
cR = zeros(out_dim,1);
Gq = zeros(out_dim,sum(1:in_ng));
Gh = zeros(out_dim,out_dim);

% Generators
for q=1:out_dim
   HessianTriangHull{q} = Interval(HessianTriangHull{q}); % Interval call just in case the result is not an interval, so mid and diam still work
    
   Qbar{q} = Z.G.'*HessianTriangHull{q}*Z.G;
   MidQbar = mid(Qbar{q});
   RadQbar = rad(Qbar{q});
   
   % cR
   cR(q) = (1/2)*trace(MidQbar);
     
   % Gq
   MidQdouble = MidQbar + MidQbar.';
   Diag = (1/4)*diag(MidQdouble).';
  
   UppTriang = zeros(1,sum(1:size(Z.G,2)-1)); counter = 0;
   for j=2:size(Z.G,2)
       for i = 1:(j-1)
           
           counter = counter + 1;
           UppTriang(counter) = MidQdouble(i,j);
       
       end
   end 
   
   Gq(q,:) = [Diag, UppTriang];
   
   % Gh
   Gh(q,q) = sum(sum(abs(RadQbar)));
    
end

% Constraints
Atildezeta = zeros(sum(1:in_nc),in_ng);
Atildecsi = zeros(sum(1:in_nc),sum(1:(in_ng-1)));
btilde = zeros(sum(1:in_nc),1);

counterROW = 0; % Row counter associated to \forall r <= s
for s=1:in_nc
    for r=1:s
        
        counterROW = counterROW + 1;       
        for i=1:in_ng
            Atildezeta(counterROW,i) = (1/2)*Z.A(r,i)*Z.A(s,i);
        end
        
        counterCOL = 0; % Column counter for Atildecsi, associated with \forall i<j
        for j=2:in_ng
            for i=1:(j-1)
                counterCOL = counterCOL + 1;
                Atildecsi(counterROW,counterCOL) = Z.A(r,i)*Z.A(s,j) + Z.A(r,j)*Z.A(s,i);               
            end
        end   
        
        btilde(counterROW) = Z.b(r)*Z.b(s) - sum(Atildezeta(counterROW,:),2); % The last sum take advantage of the summands being the same elements from Atildezeta
            
    end
end

p = Z.c - hz;
Atilde = [Atildezeta, Atildecsi, zeros(sum(1:in_nc), out_dim)];


% Remainder
Rem = CZonotope(cR, [Gq, Gh], Atilde, btilde);
if(norm(p)>1e-13) % Tolerance value for zero. If p is not zero, add the inclusion part
    L = Interval(zeros(out_dim,in_dim),zeros(out_dim,in_dim));
    for q=1:out_dim
        L(q,:) = p.'*HessianTriangHull{q};
    end    
    Rem = Rem + inclusion(CZonotope(p,2*Z.G,Z.A,Z.b),L);   
end

% -----------------------------------------------------------------------

% ---------------------- Constrained zonotope Y -------------------------

Jacobh_x = Jacobgh(:,1:nof_states);
Jacobh_v = Jacobgh(:,nof_states+1:end);
Y = (y - gh + Jacobh_x*hx + Jacobh_v*hv) - Jacobh_v*V - Rem;
%Y = CZonotope_sum(CZonotope(y - fh + Jacobh_x*hx + Jacobh_v*hv - Jacobh_v*V.c, -Jacobh_v*V.G, V.A, V.b),...
%                  CZonotope_limage(Rem,-1));
              
% -----------------------------------------------------------------------

% Generalized intersection
Xhat = intersection(X,Y,Jacobh_x);  

end



