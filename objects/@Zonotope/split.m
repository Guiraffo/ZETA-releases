function Zsplit = split(Z,gi)
%SPLIT splits a zonotope into (potentially overlapping) smaller zonotopes
%      across the directions of chosen generators
%
%   SYNTAX: Zsp = SPLIT(Z,gi)
%
%   INPUTS
%           Z: zonotope object
%          gi: array of indexes of the chosen generators
%
%   OUTPUT
%      Zsplit: the resulting list of zonotopes stored in a cell array

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

Z = Z(:);
gi = gi(:);

nof_gi = size(gi,1);

if nof_gi==1
    Zsplit{1} = Zonotope(Z.c - Z.G(:,gi)/2, Z.G); Zsplit{1}.G(:,gi) = Zsplit{1}.G(:,gi)/2;  
    Zsplit{2} = Zonotope(Z.c + Z.G(:,gi)/2, Z.G); Zsplit{2}.G(:,gi) = Zsplit{2}.G(:,gi)/2;  
    return;
else
    
    

%nof_sets = 2^nof_gi;

Zsplit_new = {Z};
for i=1:nof_gi
    
    Zsplit_prev = Zsplit_new;
    Zsplit_new = repmat({Zonotope},2^i,1);
    
    for j=1:size(Zsplit_prev,1)
        Zsplit_new{j*2-1} = Zonotope(Zsplit_prev{j}.c - Zsplit_prev{j}.G(:,gi(i))/2, Zsplit_prev{j}.G); Zsplit_new{j*2-1}.G(:,gi(i)) = Zsplit_new{j*2-1}.G(:,gi(i))/2;  
        Zsplit_new{j*2}   = Zonotope(Zsplit_prev{j}.c + Zsplit_prev{j}.G(:,gi(i))/2, Zsplit_prev{j}.G); Zsplit_new{j*2}.G(:,gi(i))   = Zsplit_new{j*2}.G(:,gi(i))/2;      
    end
    
end    

Zsplit = Zsplit_new;


end



