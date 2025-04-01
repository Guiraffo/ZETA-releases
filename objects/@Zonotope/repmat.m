function Znew = repmat(Z,N)
%REPMAT returns the cartesian product of N copies of a zonotope Z
%
%   SYNTAX: Znew = REPMAT(Z,N)
%
%   INPUTS
%           Z: zonotope
%           N: number of copies of Z 
%
%   OUTPUT
%        Znew: the cartesian product of N copies of Z

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

G_ = repmat({Z.G},N,1);
Znew = Zonotope(repmat(Z.c,N,1), blkdiag(G_{:}));

end



