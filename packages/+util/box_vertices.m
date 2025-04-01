function Vertices = box_vertices(Box,nvert)
%BOX_VERTICES returns the vertices of an interval
%
%   SYNTAX: Vertices = BOX_VERTICES(Box)
%                      BOX_VERTICES(Box,nvert)
%
%   INPUTS
%         Box: the input interval
%       nvert: (optional) number of vertices to generate
%              (algorithm will return the first nvert vertices of Box)
%              
%   OUTPUT
%    Vertices: the vertices of the input interval

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
   
space_dimension = size(Box,1);
nof_vertices = 2^space_dimension;

% Check inputs
if(nargin>1)
    if(nvert<nof_vertices)
        nof_vertices = nvert; % compute only the first nvert vertices. If nvert is bigger than 2^n, compute all the 2^n vertices.
    end
end
        

Vertices = zeros(space_dimension,nof_vertices);

c = mid(Box);
G = diag(rad(Box));
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





