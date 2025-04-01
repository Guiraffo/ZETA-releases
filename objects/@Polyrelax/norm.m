function Z = norm(X,r)
%NORM returns the norm of a polyrelax object vector
%
%   SYNTAX: Z = NORM(X)
%   SYNTAX: Z = NORM(X,r)
%
%   INPUTS
%           X: the input polyrelax object vector
%           r: 2 for Euclidean norm, or 2-norm (default)
%              1 for 1-norm
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

% Default is Euclidean norm
if(nargin==1)
    r = 2;
end

% Calls abs, or sqrt and mpower, to compute it through operator overload

X = X(:);
nof_x = length(X);

if(nof_x==1) % Scalar
    Z = abs(X);
else % Vector
    switch r
        case 2 % Euclidean norm
            Z = X(1)^2;
            for j=2:nof_x
                Z = Z + X(j)^2;
            end
            Z = sqrt(Z);
        case 1 % 1-norm
            Z = abs(X(1));
            for j=2:nof_x
                Z = Z + abs(X(j));
            end
        otherwise
            error('Invalid norm choice in Polyrelax/norm');
    end     
end

end