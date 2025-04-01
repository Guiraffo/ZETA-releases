function Znew = vertcat(varargin)
%[;] or VERTCAT returns the cartesian product of zonotopes
%
%   SYNTAX: Znew = [Z1; Z2; ... ]
%           Znew = VERTCAT(Z1,Z2,...)
%
%   INPUTS
%           Zi: the i-th zonotope, i=1,2,...
%
%   OUTPUT
%           Znew: the resulting zonotope

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

% Process inputs 
c = cell(nargin,1);
G = cell(nargin,1);
for i=1:nargin
    if(isempty(varargin{i}))
        Z = Zonotope;
    else
        Z = Zonotope(varargin{i});
    end
    c{i} = Z.c;
    G{i} = Z.G;
end
    
Znew = Zonotope(vertcat(c{:}),blkdiag(G{:}));    

end



