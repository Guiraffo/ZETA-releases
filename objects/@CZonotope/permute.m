function Znew = permute(Z,dims,gens,cons)
%PERMUTE permutes the dimensions, generators, and constraints of
%                  a constrained zonotope to a desired order
%
%   SYNTAX: Znew = PERMUTE(Z,dims,gens,cons)
%
%   INPUTS
%           Z: contrained zonotope as structure object
%        dims: array of indexes corresponding to the desired dimension
%              order (leave it empty if unchanged)
%        gens: array of indexes corresponding to the desired generator
%              order (leave it empty if unchanged)
%        cons: array of indexes corresponding to the desired constraint
%              order (leave it empty if unchanged)
%
%   OUTPUT
%        Znew: the permuted constrained zonotope

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
cons = cons(:); 

Znew = Z;
if(~isempty(dims))
    Znew = CZonotope(Znew.c(dims(:),1), Znew.G(dims(:),:), Znew.A, Znew.b);
end
if(~isempty(gens))
    Znew = CZonotope(Znew.c, Znew.G(:,gens(:)), Znew.A(:,gens(:)), Znew.b);
end
if(~isempty(cons))
    Znew = CZonotope(Znew.c, Znew.G, Znew.A(cons(:),:), Znew.b(cons(:),1));
end


end



