function V = vrep(Z)
%VREP converts a zonotope into a polytope in vertex representation
%
%   SYNTAX: V = VREP(Z)
%
%   INPUTS
%           Z: zonotope object
%
%   OUTPUT
%           V: vertices of the vertex representation

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

zonotope_dimension = size(Z.G,1);
zonotope_order = size(Z.G,2);

nof_vertices = 2^zonotope_order;
box_vertex = zeros(zonotope_order,1);
zonotope_vertices = zeros(zonotope_dimension,nof_vertices);

% Zonotope vertices by affine transformation of each vertex of the unitary box
for i=0:nof_vertices-1
    for j=0:zonotope_order-1
        if bitand(i,2^j)~=2^j
            box_vertex(j+1) = -1;
        else
            box_vertex(j+1) = 1;
        end
    end
    zonotope_vertices(:,i+1) = Z.c + Z.G*box_vertex;
end

V = zonotope_vertices;

end

