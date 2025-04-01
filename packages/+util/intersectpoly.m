function [H,k,A,b] = intersectpoly(H1,k1,A1,b1,H2,k2,A2,b2)
%INTERSECTPOLY intersects two polyhedrons in halfspace representation
%              {Hx<=k, Ax=b}
%
%   SYNTAX: INTERSECTPOLY(H1,k1,A1,b1,H2,k2,A2,b2)
%
%   INPUTS
%    H1,k1,A1,b1: halfspace representation variables of the first
%                 polyhedron
%    H2,k2,A2,b2: halfspace representation variables of the second
%                 polyhedron
%
%   OUTPUT
%        H,k,A,b: halfspace representation of the intersection

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

% Check inputs
dim1 = size(H1,2);
ins1 = size(H1,1);
dim2 = size(H2,2);
ins2 = size(H2,1);

% Empty equality variables
if(isempty(A1))
    A1 = zeros(0,dim1);
end
if(isempty(b1))
    b1 = zeros(0,1);
end
if(isempty(A2))
    A2 = zeros(0,dim2);
end
if(isempty(b2))
    b2 = zeros(0,1);
end

eqs1 = size(A1,1);
eqs2 = size(A2,1);

% Dimensions
if((dim1~=dim2)||(size(k1,1)~=ins1)||(size(k1,2)~=1)||(size(A1,2)~=dim1)||(size(b1,1)~=eqs1)||(size(b1,2)~=1)||...
                 (size(k2,1)~=ins2)||(size(k2,2)~=1)||(size(A2,2)~=dim2)||(size(b2,1)~=eqs2)||(size(b2,2)~=1))
    error('Input dimensions mismatch.')
end

H = [H1;H2];
k = [k1;k2];
A = [A1;A2];
b = [b1;b2];

end
