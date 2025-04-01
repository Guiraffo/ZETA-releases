function Znew = intersection(Z,Y,R)
%INTERSECTION returns the generalized intersection of two line
%             zonotopes
%
%   SYNTAX: Znew = INTERSECTION(Z,Y)
%           Znew = INTERSECTION(Z,Y,R)
%
%   INPUTS
%              Z: the first line zonotope as a structure
%              Y: the second line zonotope as a structure
%              R: linear mapping matrix (generalized intersection, default
%                 value is the identity matrix)
%
%   OUTPUT
%           Znew: the resulting line zonotope

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
Z = LZonotope(Z);
Y = LZonotope(Y);
setdim = Z.dim;
if(nargin<2)
    error('Not enough input arguments.');
elseif(nargin==2)
    if(setdim~=Y.dim)
        error('Wrong input dimensions in LZonotope/intersection.')
    end    
    R = eye(setdim);
elseif(nargin==3)
    if(~isnumeric(R))
        error('Invalid inputs in LZonotope/intersection.')
    elseif((size(R,2)~=setdim)||(size(R,1)~=Y.dim))
        error('Wrong input dimensions in LZonotope/intersection.')
    end
else
    error('Invalid inputs in LZonotope/intersection.');
end
    

% Intersection of two line zonotopes
       
newc = Z.c;
newG = [Z.G, zeros(size(Z.G,1),size(Y.G,2))];
newM = [Z.M, zeros(size(Z.M,1),size(Y.M,2))];    
newS = [  Z.S,                              zeros(size(Z.S,1),size(Y.S,2));
          zeros(size(Y.S,1), size(Z.S,2)),                             Y.S;
                                    R*Z.M,                            -Y.M];
newA = [  Z.A,                              zeros(size(Z.A,1),size(Y.A,2));
          zeros(size(Y.A,1), size(Z.A,2)),                             Y.A;
                                    R*Z.G,                            -Y.G];
newb = [Z.b; 
        Y.b;
        Y.c - R*Z.c];

Znew = LZonotope(newc,newG,newM,newS,newA,newb);
     

end
