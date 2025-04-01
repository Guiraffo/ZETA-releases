function result = isempty(Z)
%ISEMPTY verifies if a constrained zonotope is empty
%
%   SYNTAX: result = ISEMPTY(Z)
%
%   INPUTS
%           Z: constrained zonotope object
%
%   OUTPUT
%      result: 0 if the constrained zonotope is not empty
%              1 if the constrained zonotope is empty

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

% Input check
if((nargin<1)||isempty(Z.c))
    result = 1;
    return;
end

nof_generators = size(Z.G,2);
nof_constraints = size(Z.A,1);


f = [1; zeros(nof_generators,1)];
A = [                 -1,  zeros(1,nof_generators);
     -ones(nof_generators,1),  eye(nof_generators);
     -ones(nof_generators,1), -eye(nof_generators)];
b = [0; zeros(2*nof_generators,1)];
Aeq = [zeros(nof_constraints,1), Z.A];
beq = Z.b;

OPTIONS.Display = 'off';

[x,~,exitflag] = optim.solvelp(f,A,b,Aeq,beq,[],[],OPTIONS);


if(ischar(exitflag)) % Gurobi
    if(strcmp(exitflag,'OPTIMAL')) % Success    norm_csi = x(1);
        norm_csi = x(1);
        if norm_csi <= 1
            result = 0;
        else
            result = 1;
        end
    else % Not success (treated as somehow unfeasible and the point z not belonging to the constrained zonotope Z)
        result = 1; 
    end
elseif(exitflag==1) % Linprog success
    norm_csi = x(1);
    if norm_csi <= 1
        result = 0;
    else
        result = 1;
    end
else % Not success (treated as somehow unfeasible and the point z not belonging to the constrained zonotope Z)
    result = 1;
end


end

