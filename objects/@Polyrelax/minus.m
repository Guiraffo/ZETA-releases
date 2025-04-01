function Z = minus(X,Y)
%+ or PLUS returns the difference of two Polyrelax objects
%
%   SYNTAX: Z = X-Y
%           Z = MINUS(X,Y)
%
%   INPUTS
%           Z: the first Polyrelax object
%           W: the second Polyrelax object
%
%   OUTPUT
%           Z: the resulting Polyrelax object

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

X = Polyrelax(X);
Y = Polyrelax(Y);

Z = Polyrelax(X.x-Y.x);

% Adds an equality constraint to Polyrelax Hrep
Aadd = zeros(1,Z.i);
badd = zeros(1,1);

if(X.i==-1) % If X is a constant
    
    % [ ... 1 ... 1]z = x
    Aadd(Y.i) = 1;
    Aadd(Z.i) = 1;    
    badd = X.x;
    
elseif(Y.i==-1) % If Y is a constant
    
    % [ ... -1 ... 1]z = -y
    Aadd(X.i) = -1;
    Aadd(Z.i) = 1;    
    badd = -Y.x;    
    
else % If neither X or Y is a constant

    % [... -1 ... 1 ... 1]z = 0
    Aadd(X.i) = -1;
    Aadd(Y.i) = 1;
    Aadd(Z.i) = 1;

end

Polyrelax.addHrepArow(Aadd,badd);
Polyrelax.addelimind(Z.i);

% Revision 05-03-2024: added proper treatment for constants
% First version: 05-02-2024
    
end