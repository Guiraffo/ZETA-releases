function Znew = convexhull(varargin)
%CONVEXHULL computes the convex hull of several zonotopes using the method 
%           in Raghuraman and Koeln (2022), iteratively. The result is a
%           constrained zonotope.
%
%   SYNTAX: Znew = CONVEXHULL(Z1,Z2,...)
%
%   INPUTS
%             Zi: the i-th zonotope
%
%   OUTPUT
%           Znew: the constrained zonotope enclosure

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

if(nargin<1)
    error('Not enough input arguments.');
end

Z = varargin;
nof_sets = length(Z);

if(nof_sets==1)
    Znew = Z;
    return;
end

% Work with structs in this step to mitigate overhead
Z{1} = CZonotope(Z{1});
Znew_.c = Z{1}.c;
Znew_.G = Z{1}.G;
Znew_.A = Z{1}.A;
Znew_.b = Z{1}.b;
for i=2:nof_sets
    Znew_ = convexhull_Raghuraman(Znew_,CZonotope(Z{i}));
end

% Resulting constrained zonotope
Znew = CZonotope(Znew_.c,Znew_.G,Znew_.A,Znew_.b);

end

function Znew = convexhull_Raghuraman(Z1,Z2)
% Computes the convex hull of two constrained zonotopes using Theorem 5 in
% 'Set operations and order reductions for constrained zonotopes', 
%  Raghuraman and Koeln (2022)

ndim1 = size(Z1.c,1);
ndim2 = size(Z2.c,1);
if(ndim1~=ndim2)
    error('Wrong dimensions in czonotope/convexhull.');
end

ng1 = size(Z1.A,2);
ng2 = size(Z2.A,2);
nc1 = size(Z1.A,1);
nc2 = size(Z2.A,1);


A31 = [eye(ng1); -eye(ng1); zeros(ng2,ng1); zeros(ng2,ng1)];
A32 = [zeros(ng1,ng2); zeros(ng1,ng2); eye(ng2); -eye(ng2)];
A30 = 0.5*[-ones(2*ng1,1); ones(2*ng2,1)];

Znew.c = 0.5*(Z1.c + Z2.c);
Znew.G = [Z1.G, Z2.G, 0.5*(Z1.c - Z2.c), zeros(ndim1,2*(ng1+ng2))];
Znew.A = [blkdiag(Z1.A,Z2.A), 0.5*[-Z1.b; Z2.b], zeros(nc1+nc2,2*(ng1+ng2));
                    A31, A32,               A30,           eye(2*(ng1+ng2))];
Znew.b = 0.5*[ Z1.b; Z2.b; -ones(2*(ng1+ng2),1)];                
         

end






