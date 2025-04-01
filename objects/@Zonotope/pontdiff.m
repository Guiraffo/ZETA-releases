function [Hnew,knew] = pontdiff(Z,H,k)
%PONTDIFF returns the Pontryagin difference P - Z where P = {Hx<=k} is a
%                   convex polytope in H-rep and Z is a zonotope
%
%   SYNTAX: Znew = PONTDIFF(Z,H,k)
%
%   INPUTS
%         H,k: polyhedron in H-rep 
%           Z: zonotope in G-rep (zonotope object)
%
%   OUTPUT
%   Hnew,knew: resulting convex polytope in H-rep

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


knew = zeros(size(k,1),1);

for i=1:size(k,1)
    
    knew(i) = k(i) - supportfun(Z,H(i,:).');
    
end

Hnew = H;      

end



