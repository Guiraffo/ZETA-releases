function Zsum = plus(Z,W)
%+ or PLUS returns the Minkowski sum of two line zonotopes
%
%   SYNTAX: Zsum = Z+W
%           Zsum = PLUS(Z,W)
%
%   INPUTS
%           Z: the first line zonotope
%           W: the second line zonotope
%
%   OUTPUT
%           Zsum: the resulting line zonotope

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

Z = LZonotope(Z);
W = LZonotope(W);

if(Z.dim~=W.dim)
    error('Inputs in + have incompatible dimensions.');
else
    Zsum = LZonotope(Z.c + W.c, [Z.G, W.G], [Z.M, W.M], blkdiag(Z.S,W.S), blkdiag(Z.A,W.A), [Z.b;W.b]);
end      

end



