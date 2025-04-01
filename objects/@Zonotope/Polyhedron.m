function Poly = Polyhedron(Z,method)
%POLYHEDRON converts a zonotope into a MPT Polyhedron object
%
%   SYNTAX: Poly = POLYHEDRON(Z)
%           Poly = POLYHEDRON(Z,method)
%
%   INPUTS
%           Z: zonotope object
%      method: 'Hrep': converts Z to half-space representation 
%              'Vrep': converts Z to vertex representation 
%              (default method is 'Hrep')
%
%   OUTPUT
%        Poly: Polyhedron object corresponding to the zonotope

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
elseif(nargin==1)
    method = 'Hrep';
end

switch method
    case 'Hrep'
        [H,k] = hrep(Z);
        Poly = Polyhedron(H,k);     
    case 'Vrep'        
        V = vrep(Z);
        Poly = Polyhedron(V.');
    otherwise       
        error('Invalid method in Zonotope/Polyhedron.')
end

end

