function Znew = reduction(Z,ng_max,nc_max)
%REDUCTION performs complexity reduction of a line zonotope
%
%   SYNTAX: Znew = REDUCTION(Z,ng_max,nc_max)
%
%   INPUTS
%           Z: line zonotope object
%      ng_max: desired number of generators (empty if no generator reduction)
%      nc_max: desired number of constraints (empty if no constraint elimination)
%
%   OUTPUT
%        Znew: the new line zonotope


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

% Performs line elimination, then constraint elimination, then generator reduction

if(nargin<3)
    error('Not enough input arguments.');
end

% Eliminate all possible lines
while(sum(sum(abs(Z.S))) > 0)
    Z = reduction_lin(Z);
end

% Constraint elimination
if(~isempty(nc_max))
    Z = reduction_con(Z,nc_max);
end

% Generator reduction
if(~isempty(ng_max))
    Z = reduction_gen(Z,ng_max);
end
Znew = Z;


end

function Znew = reduction_con(Z,nc_max)
%REDUCTION_CON performs constraint reduction of a line zonotope
%
%   SYNTAX: Znew = REDUCTION_CON(Z,nc_max)
%
%   INPUTS
%           Z: line zonotope object
%      nc_max: desired number of constraints
%
%   OUTPUT
%        Znew: the new line zonotope

if(~(isempty(Z.S) || sum(sum(abs(Z.S)))==0)) % If there is any cross constraint
    disp('Cross constraint detected in LZonotope/reduction_con. Constraint elimination could not proceed.')
    Znew = Z;
    return;
end

% Constraint elimination of the bounded part
Zaux = CZonotope(Z.c, Z.G, Z.A, Z.b);
Zaux = reduction(Zaux,[],nc_max);

Znew = LZonotope(Zaux.c, Zaux.G, Z.M, zeros(size(Zaux.A,1),size(Z.M,2)), Zaux.A, Zaux.b);

end

function Znew = reduction_gen(Z,ng_max)
%REDUCTION_GEN performs generator reduction of a line zonotope. Possible
%              only for line zonotopes without cross constraints between
%              lines and generators.
%
%   SYNTAX: Znew = GREDUCTION_GEN(Z,ng_max)
%
%   INPUTS
%           Z: line zonotope object
%      ng_max: desired number of generators
%
%   OUTPUT
%        Znew: the new line zonotope


if(isempty(Z.S) || sum(sum(abs(Z.S)))) % Verifies if there is any cross constraint
    
    Zaux = CZonotope(Z.c, Z.G, Z.A, Z.b);
    Zaux = reduction(Zaux,ng_max,[]); % The unbounded part should not affect the generator reduction of the bounded part
    Znew = LZonotope(Zaux.c, Zaux.G, Z.M, Z.S, Zaux.A, Zaux.b);    
       
else   
    disp('Cross constraint detected in LZonotope/reduction_gen. Generator reduction cannot proceed.');
    Znew = Z;
end

end

function Z = reduction_lin(Z)
%REDUCTION_LIN removes one line from a line zonotope by solving one
%              constraint for the respective line
%
%   SYNTAX: Znew = REDUCTION_LIN(Z)
%
%   INPUTS
%           Z: line zonotope object
%      nc_max: desired number of constraints
%
%   OUTPUT
%        Znew: the new line zonotope



[~,maxindex] = max(abs(Z.S(:,1)));
smax = Z.S(maxindex,1);

tol = max(size(Z.S,1),size(Z.S,2))*eps('double')*norm(Z.S,'inf'); % Tolerance value

n = size(Z.c,1);
nc = size(Z.A,1);

if(abs(smax) >= tol)
    Z.A(maxindex,:) = Z.A(maxindex,:)/smax;
    Z.S(maxindex,:) = Z.S(maxindex,:)/smax;
    Z.b(maxindex,:) = Z.b(maxindex,:)/smax;
    
    xi_elim = 1;
    pivot_row = maxindex;
    
    for row=1:pivot_row-1
        lambda = Z.S(row,xi_elim);
        Z.b(row  ) = Z.b(row  )-lambda*Z.b(pivot_row  ); 
        Z.A(row,:) = Z.A(row,:)-lambda*Z.A(pivot_row,:); 
        Z.S(row,:) = Z.S(row,:)-lambda*Z.S(pivot_row,:); 
    end

    for row=pivot_row+1:nc
        lambda = Z.S(row,xi_elim);
        Z.b(row  ) = Z.b(row  )-lambda*Z.b(pivot_row  ); 
        Z.A(row,:) = Z.A(row,:)-lambda*Z.A(pivot_row,:); 
        Z.S(row,:) = Z.S(row,:)-lambda*Z.S(pivot_row,:); 
    end

    for row=1:n
        lambda = Z.M(row,xi_elim);
        Z.c(row  ) = Z.c(row  )+lambda*Z.b(pivot_row  ); 
        Z.G(row,:) = Z.G(row,:)-lambda*Z.A(pivot_row,:); 
        Z.M(row,:) = Z.M(row,:)-lambda*Z.S(pivot_row,:); 
    end      
   
    Z.M(:,xi_elim) = [];
    Z.S(:,xi_elim) = [];
    Z.A(pivot_row,:) = [];
    Z.S(pivot_row,:) = [];
    Z.b(pivot_row,:) = []; 
    
else
    
    error('Error in LZonotope/reduction_lin');
    
end


end
