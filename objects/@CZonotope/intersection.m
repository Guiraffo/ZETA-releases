function Znew = intersection(Z,varargin)
%INTERSECTION returns the generalized intersection of two constrained
%             zonotopes, or the intersection of a constrained zonotope
%             and a convex polytope in halfspace represesntation
%
%   SYNTAX: Znew = INTERSECTION(Z,Y)
%           Znew = INTERSECTION(Z,Y,R)
%           Znew = INTERSECTION(Z,Hineq,kineq)
%           Znew = INTERSECTION(Z,Hineq,kineq,Aeq,beq)
%
%   INPUTS
%              Z: the first contrained zonotope as a structure
%              Y: the second contrained zonotope as a structure
%              R: linear mapping matrix (generalized intersection, default
%                 value is the identity matrix)
%    Hineq,kineq: inequality constraints of the convex polytope
%        Aeq,beq: equality constraints of the convex polytope
%
%   OUTPUT
%           Znew: the resulting constrained zonotope

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

% Input check

intCZCZ = 0; % Flag for CZ-CZ intersection. 0 means CZ-Poly intersection
setdim = Z.dim;
if(nargin<2)
    error('Not enough input arguments.');
elseif(nargin==2)
    Y = varargin{1};
    if(isa(Y,'CZonotope'))
        intCZCZ = 1;        
        if(setdim~=Y.dim)
            error('Wrong input dimensions in CZonotope/intersection.')
        end    
        R = eye(setdim);
    else
        error('Invalid inputs in CZonotope/intersection.')
    end
elseif(nargin==3)
    Y = varargin{1};
    if(isa(Y,'CZonotope'))
        intCZCZ = 1;
        R = varargin{2};
        if(~isnumeric(R))
            error('Invalid inputs in CZonotope/intersection.')
        elseif((size(R,2)~=setdim)||(size(R,1)~=Y.dim))
            error('Wrong input dimensions in CZonotope/intersection.')
        end
    elseif(isnumeric(Y)&&isnumeric(varargin{2}))
        Hineq = Y;
        kineq = varargin{2};
        Aeq = zeros(0,setdim);
        beq = zeros(0,1);        
        if((size(Hineq,2)~=setdim)||(size(Hineq,1)~=size(kineq,1))||(size(kineq,2)~=1))
            error('Wrong input dimensions in CZonotope/intersection.')
        end
    else
        error('Invalid inputs in CZonotope/intersection.')
    end
elseif(nargin==5)
    Hineq = varargin{1};
    kineq = varargin{2};
    Aeq = varargin{3};
    beq = varargin{4};
    if(~(isnumeric(Hineq)&&isnumeric(kineq)&&isnumeric(Aeq)&&isnumeric(beq)))
        error('Invalid inputs in CZonotope/intersection.');
    elseif((size(Hineq,2)~=setdim)||(size(Hineq,1)~=size(kineq,1))||(size(kineq,2)~=1)||(size(Aeq,2)~=setdim)||(size(Aeq,1)~=size(beq,1))||(size(beq,2)~=1))
        error('Wrong input dimensions in CZonotope/intersection.')
    end
else
    error('Invalid inputs in CZonotope/intersection.');
end
    

if(intCZCZ) % Intersection of two constrained zonotopes
       
    newc = Z.c;
    newG = [Z.G, zeros(size(Z.G,1),size(Y.G,2))];
    newA = [  Z.A,                              zeros(size(Z.A,1),size(Y.A,2));
              zeros(size(Y.A,1), size(Z.A,2)),                             Y.A;
                                        R*Z.G,                            -Y.G];
    newb = [Z.b; 
            Y.b;
            Y.c - R*Z.c];
      
    Znew = CZonotope(newc,newG,newA,newb);
      
else % Intersection of a constrained zonotope and a convex polytope
    
    Znew = intersection_hrep(Z,Hineq,kineq,Aeq,beq);
    
end

end


function Znew = intersection_hrep(Z,H,k,Aeq,beq)
%INTERSECTION_HREP returns the intersection of a constrained zonotope and a
%                  polytope in H-rep (the latter is not necessarily bounded)
%
%   SYNTAX: Znew = INTERSECTION_HREP(Z,P)
%
%   INPUTS
%              Z: constrained zonotope
%            H,k: inequality constraints of the convex polytope
%        Aeq,beq: equality constraints of the convex polytope
%
%   OUTPUT
%           Znew: the resulting constrained zonotope


% Lower bound of the interval hull of H*Z
sigma = H*Z.c - sum(abs(H*Z.G),2);


newc = Z.c;
newG = [Z.G, zeros(size(Z.G,1),size(H,1))];
newA = [Z.A, zeros(size(Z.A,1),size(H,1));
        H*Z.G, diag(sigma-k)/2;
        Aeq*Z.G, zeros(size(Aeq,1),size(k,1))];
newb = [Z.b;
        (sigma+k)/2 - H*Z.c;
        beq - Aeq*Z.c];
      
Znew = CZonotope(newc,newG,newA,newb);

end







