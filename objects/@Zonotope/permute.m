function Znew = permute(Z,dims,gens)
%PERMUTE permutes the dimensions and generators of a zonotope to a desired
%        order
%
%   SYNTAX: Znew = PERMUTE(Z,dims,gens)
%
%   INPUTS
%           Z: zonotope object
%        dims: array of indexes corresponding to the desired dimension
%              order (leave it empty if unchanged)
%        gens: array of indexes corresponding to the desired generator
%              order (leave it empty if unchanged)
%
%   OUTPUT
%        Znew: the permuted zonotope

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

dims = dims(:);
gens = gens(:);

Znew = Z;
if(~isempty(dims))
    Znew = Zonotope(Znew.c(dims(:),1), Znew.G(dims(:),:));
end
if(~isempty(gens))
    Znew = Zonotope(Znew.c, Znew.G(:,gens(:)));
end


% First version: 12-05-2020

end



