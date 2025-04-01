function P = partopebound(Z,bound_method,red_method)
%PARTOPEBOUND returns a parallelotope enclosing a zonotope
%
%    SYNTAX: P = PARTOPEBOUND(Z)
%            P = PARTOPEBOUND(Z,bound_method)
%            P = PARTOPEBOUND(Z,bound_method,red_method)
%
%    INPUTS
%            Z: zonotope object
%              
%    OUTPUT
%            P: parallelotope enclosing the zonotope Z (P is a zonotope object)
% bound_method: - 'SVD': uses Lemma 3 in Alamo et al (2008)
%               - 'Reduction': uses generator reduction (according to red_method)
%               - 'Interval': uses the interval hull of Z
%              Default bound_method is 'SVD'.
%   red_method:- 'Combastel': uses Method 1 in Yang and Scott (2018)
%              - 'Girard': uses Method 2 in Yang and Scott (2018) 
%              - 'Chisci': uses Method 4 in Yang and Scott (2018)
%              Default red_method is 'Chisci'.

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

switch nargin
    case 0
        error('Not enough input arguments.');
    case 1
        bound_method = 'SVD';
    case 2
        red_method = 'Chisci';
end


switch bound_method
    case 'SVD'

        [U,S,V] = svd(Z.G);
        P = Zonotope(Z.c,U*diag(sum(abs(S*V.'),2)));

    case 'Reduction'
        
        P = reduction(Z,Z.dim,red_method);
        
    case 'Interval'
        
        P = Zonotope(intervalhull(Z));

    otherwise       
        error('Invalid bounding method in partopebound.');
end

end



