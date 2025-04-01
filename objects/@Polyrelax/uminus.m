function Z = uminus(X)
% - Unary minus returns the negative of a polyrelax object
%
%   SYNTAX: Z = -X
%           Z = UMINUS(X)
%
%   INPUTS
%           X: the polyrelax object
%
%   OUTPUT
%           Z: the resulting polyrelax object

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

Z = Polyrelax(-X.x);


% Adds an equality constraint to polyrelax Hrep
% [... 1 ... 1]z = 0
Aadd = zeros(1,Z.i);
badd = zeros(1,1);

Aadd(X.i) = 1;
Aadd(Z.i) = 1;

Polyrelax.addHrepArow(Aadd,badd);
Polyrelax.addelimind(Z.i);

% First version: 05-02-2024
    
end