function Vertices = partope_vertices(Partope)
%PARTOPE_VERTICES returns the vertices of a parallelotope
%
%   SYNTAX: Vertices = PARTOPE_VERTICES(Partope)
%
%   INPUTS
%     Partope: the input paralelotope (zonotope object)
%              
%   OUTPUT
%    Vertices: the vertices of Partope

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
   
space_dimension = size(Partope.c,1);
nof_vertices = 2^space_dimension;
Vertices = zeros(space_dimension,nof_vertices);

c = Partope.c;
G = Partope.G;
unitbox_vertex = zeros(space_dimension,1);

% Box vertices by affine transformation of each vertex of the unitary box
for i=0:nof_vertices-1
    
    for j=0:space_dimension-1
        if bitand(i,2^j)~=2^j
            unitbox_vertex(j+1) = -1;
        else
            unitbox_vertex(j+1) = 1;
        end
    end
    Vertices(:,i+1) = c + G*unitbox_vertex;
end        

end





