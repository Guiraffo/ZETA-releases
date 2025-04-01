function Znew = propagate(Z,fname,OPTIONS)
%PROPAGATE propagates a zonotope through a nonlinear function using a 
%          chosen approximation method
%
%   SYNTAX: Znew = propagate(Z,fname)
%           Znew = propagate(Z,fname,OPTIONS)
%   INPUTS
%           Z: zonotope to be propagated
%       fname: function name of the nonlinear function
%     OPTIONS: (optional) OPTIONS structure claimed using util.propagopts
%              (default value is util.propagopts('MV'))
%
%   OUTPUT
%         Znew: zonotope enclosing the propagated Z through f

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

if(nargin<2)
    error('Not enough input arguments.');
elseif(nargin==2)
    OPTIONS = util.propagopts('ZMV');
end

method = OPTIONS.method;
f = util.fhandle(fname);

switch method
    case 'ZMV'
        Znew = propagate_MV(Z,f);
    case 'ZFO'
        Znew = propagate_FO(Z,f);
    case 'ZDC'
        Znew = propagate_DC(Z,f,OPTIONS);
    otherwise
        error('Invalid propagation method in Zonotope/propagate.');
end

end

%% Mean value extension
function Znew = propagate_MV(Z,f)
% Mean value extension based on Theorem 4 from "Guaranteed State Estimation by Zonotopes" (Alamo et al, 2005)

% Interval hull and Jacobian
Zhull = intervalhull(Z);
Jacobf = f.Jacob(Zhull);

% Interval matrix M
M = Jacobf*Z.G;

% Zonotope inclusion and propagated set
Zinc = Zonotope.inclusion(zeros(size(M,1),1), M);
Znew = f.eval(Z.c) + Zinc;

end

%% First-order Taylor extension
function Znew = propagate_FO(Z,f)
% First-order Taylor extension based on Combastel (2005) 

% Interval hull
Zhull = intervalhull(Z);

% Evaluate the nonlinear function, Jacobian, and triangular Hessians

fc = f.eval(Z.c);
Jacobfc = f.Jacob(Z.c);
HessianTriangHull = f.HessianTriang(Zhull);   
    

% Initialize auxiliary variables
in_dim = Z.dim;
in_ng = Z.ng;
out_dim = size(fc,1);
Qbar = cell(out_dim,1);
cR = zeros(out_dim,1);
Gq = zeros(out_dim,sum(1:in_ng));
Gh = zeros(out_dim,out_dim);

for q=1:out_dim
    
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

Znew = Zonotope(fc + cR, [Jacobfc*Z.G, Gq, Gh]); 

end

%% DC programming principles
function Znew = propagate_DC(Z,f,OPTIONS)
% Propagation based on DC programming principles (Alamo et al., 2008)

if(OPTIONS.partopeimprove)
    [Zbar,P_vert,hz,lambda] = dcprog.propagate(Z,f,OPTIONS);
    Znew = dcprog.partopeimprove(Zbar,f,P_vert,hz,OPTIONS,lambda);
else
    Znew = dcprog.propagate(Z,f,OPTIONS);
end

end
