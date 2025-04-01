function [Xhat,Xbar,cputimes] = zon_firstorder_estimator(dtsys,Xhat_prev,W,V,u,y,OPTIONS)
%ZON_FIRSTORDER_ESTIMATOR implements a state estimator for nonlinear
%                         discrete-time systems using zonotopes based on
%                         the first-order Taylor extension

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
Xhat = reduction(Xhat,OPTIONS.ngmax);
cputimes.redu = toc;

% Total computational time
cputimes.total = cputimes.pred + cputimes.updt + cputimes.redu;

end

function Xbar = prediction(f,X,W,u,Ts)
% Implementation of the first order taylor extension in Combastel (2005)
% for state estimation

% Interval hulls
Xhull = intervalhull(X);
if(isempty(W))
    W = Zonotope;
    Whull = [];
else
    Whull = intervalhull(W);
end

Z = [X;W];

% Evaluate the nonlinear function, Jacobian, and triangular Hessians

fc = f.eval(X.c,W.c,u,Ts);
Jacobfc = f.Jacob(X.c,W.c,u,Ts);
HessianTriangHull = f.HessianTriang(Xhull,Whull,u,Ts);

% Initialize auxiliary variables
in_dim = Z.dim;
in_ng = Z.ng;
out_dim = size(fc,1);
Qbar = cell(out_dim,1);
cR = zeros(out_dim,1);
Gq = zeros(out_dim,sum(1:in_ng));
Gh = zeros(out_dim,out_dim);

for q=1:out_dim
   HessianTriangHull{q} = interval(HessianTriangHull{q}); % Interval call just in case the result is not an interval, so mid and diam still work  
    
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

Xbar = Zonotope(fc + cR, [Jacobfc*Z.G, Gq, Gh]); 

end

function Xhat = update(g,Xbar,V,y)

Xtilde = Xbar;
if(isempty(V))
    V = Zonotope;
    Vhull = [];
else
    Vhull = intervalhull(V);
end

% For loop

for i=1:length(y)
    
    % Compute the strip of consistent states
    
    % Get interval Jacobian and the strip parameter p
    Xhull = intervalhull(Xtilde);
    Jacobg = g.Jacob(Xhull,Vhull);
    Jacobg_x = Jacobg(:,1:size(Xhull,1));

    p = mid(Jacobg_x(i,:).');

    % Computation of the strip parameters d and sigma

    if(sum(diam(Jacobg(i,:)))==0) % If the Jacobian is not an interval, it means this measurement is linear

        Dvi = mid(Jacobg(i,size(Xhull,1)+1:end));
        DvVi = Dvi.'*V;
        d = y(i) - DvVi.c(i);
        sigma = sum(rad(intervalhull(DvVi)));

    else % Otherwise, use the general formula in Alamo et al (2005)

        outputhull = model.g(Xhull,Vhull);
        stripinterval = p.'*Xhull - outputhull(i);
        d = y(i) + mid(stripinterval);
        sigma = rad(stripinterval);

    end
    
    % Intersect with the strip of consistent states
    Xtilde = intersection(Xtilde,strip(p, d, sigma));
    
end
    
Xhat = Xtilde;

end

