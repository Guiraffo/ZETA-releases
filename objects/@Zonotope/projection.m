function Znew = projection(Z,axes)
%_PROJECTION returns the projection of a zonotope along a given set of axes
%
%   SYNTAX: Znew = PROJECTION(Z,axes)
%
%   INPUTS
%           Z: zonotope object
%        axis: array of indexes associated to the axes of interest
%
%   OUTPUT
%        Znew: the projected zonotope

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

dim = Z.dim;
dim_proj = length(axes);
proj = zeros(dim_proj, dim);

for j=1:dim_proj
    proj(j,axes(j)) = 1;
end

Znew = proj*Z;

end



