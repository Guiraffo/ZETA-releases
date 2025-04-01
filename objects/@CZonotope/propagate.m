function Znew = propagate(Z,fname,OPTIONS)
%PROPAGATE propagates a constrained zonotope through a nonlinear function
%          using a chosen approximation method
%
%   SYNTAX: Znew = propagate(Z,fname)
%           Znew = propagate(Z,fname,OPTIONS)
%   INPUTS
%           Z: constrained zonotope to be propagated
%       fname: function name of the nonlinear function
%     OPTIONS: (optional) OPTIONS structure claimed using util.propagopts
%              (default value is util.propagopts('CZMV'))
%
%   OUTPUT
%        Znew: constrained zonotope enclosing the propagated Z through f

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
    OPTIONS = util.propagopts('CZMV');
end

method = OPTIONS.method;
f = util.fhandle(fname);

switch method
    case 'CZMV'
        Znew = propagate_MV(Z,f);
    case 'CZFO'
        Znew = propagate_FO(Z,f);
    case 'CZDC'
        Znew = propagate_DC(Z,f,OPTIONS);
    case 'CZPR'
        Znew = propagate_PR(Z,f,OPTIONS);        
    otherwise
        error('Invalid OPTIONS in CZonotope/propagate.');
end

end

%% Mean value extension
function Znew = propagate_MV(Z,f)
% Mean value extension using constrained zonotopes (Rego et al., 2021)

% Interval hull and Jacobian
Zhull = intervalhull(Z);
f_Jacob = f.Jacob(Zhull);

% Get approximation point h
[h,~,Zbar] = chooseh(Z,Zhull,f_Jacob,'mindiamm');

% Zonotope inclusion and propagated set
Zinc = inclusion(CZonotope(Z.c-h,Z.G,Z.A,Z.b),f_Jacob,Zbar);
Znew = f.eval(h) + Zinc;

end

%% First-order Taylor extension
function Znew = propagate_FO(Z,f)
% First-order Taylor extension (Rego et al., 2021)

Zhull = intervalhull(Z);

% Choose h
h = chooseh(Z,Zhull,[],'distcenter');

% Evaluate the nonlinear function, Jacobian, and triangular Hessians

fh = f.eval(h);
Jacobfh = f.Jacob(h);
HessianTriangHull = f.HessianTriang(Zhull);   
    

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

p = Z.c - h;
Atilde = [Atildezeta, Atildecsi, zeros(sum(1:in_nc), out_dim)];
Znew = CZonotope(fh + Jacobfh*p + cR, [Jacobfh*Z.G, Gq, Gh],...
                  [Z.A, zeros(in_nc,size(Gq,2)+size(Gh,2)); zeros(size(Atilde,1), in_ng), Atilde],...
                  [Z.b; btilde]); 

if(norm(p)>1e-13) % Tolerance value for zero. If p is not zero, add the inclusion part

    L = Interval(zeros(out_dim,in_dim),zeros(out_dim,in_dim));
    
    for q=1:out_dim
        L(q,:) = p.'*HessianTriangHull{q};
    end
    
    Znew = Znew + inclusion(CZonotope(p,2*Z.G,Z.A,Z.b),L);
    
end

% % -----------------------------------------------------------------------

end

%% DC programming principles
function Znew = propagate_DC(Z,f,OPTIONS)
% Propagation based on DC programming principles (de Paula et al., 2024)

Znew = dcprog.propagate(Z,f,OPTIONS);

end

%% Polyhedron relaxations
function Snew = propagate_PR(S,f,OPTIONS)
% Propagation based on polyhedron relaxations (Rego et al., 2024)

% Clear and initialize the static variables hrep, z, and elimind
Polyrelax.clear;

% Initialize Polyrelax objects for the input
Shull = intervalhull(S);
nof_s = S.dim;
nof_inps = nof_s;
Spolyrel = Polyrelax.empty(0,1);
for j=1:nof_s
    Spolyrel(j,1) = Polyrelax(Shull(j));
end

% Propagate the Polyrelax objects through the nonlinear function f
Zpolyrel = f.eval(Spolyrel);

% Get the indexes of the new X (or old X if not prediction)
projindexes = [Zpolyrel.i];

% Partially solve A*z=b in the H-rep and update projindexes
if(OPTIONS.solvequalities)
    [Gh,ch,inputind] = Polyrelax.solveAbTriangPartial(projindexes,nof_inps);
else
    inputind = 1:nof_inps;
end

% Get resulting Hrep and Z, plus projection indexes of the image of f
thisHrep = Polyrelax.Hrep;
thisZ = Polyrelax.Z;
nof_inps = length(inputind);

% Augmented constrained zonotope S x Ztil
Saug = [projection(S,inputind); CZonotope(thisZ(nof_inps+1:end))];

% Intersection of (S x Ztil) with the polyhedral relaxation
Saug_n_P = intersection(Saug,thisHrep.H,thisHrep.k,thisHrep.A,thisHrep.b);

% Projects the result into the image of f
if(OPTIONS.solvequalities)
    Snew = Gh*Saug_n_P + ch;    
else
    Snew = projection(Saug_n_P,projindexes);
end


end
