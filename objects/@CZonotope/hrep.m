function [H,k] = hrep(Z)
%HREP converts a constrained zonotope into halfspace representation
%
%   SYNTAX: [H,k] = HREP(Z)
%
%   INPUTS
%           Z: constrained zonotope object
%
%   OUTPUT
%       [H,k]: constrained zonotope in halfspace representation

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

if(nargin<1)
    error('Not enough input arguments.');
end

    
[Zlift,dimension] = lift(Z);
[H_,k] = hrep(Zlift);
H = H_(:,1:dimension);
          
end

