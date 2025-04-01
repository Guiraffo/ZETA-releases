function result = isempty(Z)
%ISEMPTY verifies if a line zonotope is empty
%
%   SYNTAX: result = ISEMPTY(Z)
%
%   INPUTS
%           Z: line zonotope object
%
%   OUTPUT
%      result: 0 if the line zonotope is not empty
%              1 if the line zonotope is empty

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

nof_generators = size(Z.G,2);
nof_linees = size(Z.M,2);
nof_constraints = size(Z.A,1);

f = [1; zeros(nof_generators,1); zeros(nof_linees,1)];
A = [-ones(nof_generators,1),  eye(nof_generators), zeros(nof_generators,nof_linees);
     -ones(nof_generators,1), -eye(nof_generators), zeros(nof_generators,nof_linees)];
b = zeros(2*nof_generators,1);
Aeq = [zeros(nof_constraints,1), Z.A, Z.S];
beq = Z.b;
       
       

% OPTIONS
persistent OPTIONS; % Helps computational time
if(isempty(OPTIONS))
    OPTIONS.Display = 'off';
end

       
[x,~,exitflag] = optim.solvelp(f,A,b,Aeq,beq,[],[],OPTIONS);

if(ischar(exitflag)) % Gurobi
    if(strcmp(exitflag,'OPTIMAL')) % Success    norm_csi = x(1);
        norm_csi = x(1);
        if norm_csi <= 1
            result = 0;
        else
            result = 1;
        end
    else % Not success (treated as somehow unfeasible and the point z not belonging to the line zonotope Z)
        result = 1; 
    end
elseif(exitflag==1) % linprog success
    norm_csi = x(1);
    if norm_csi <= 1
        result = 0;
    else
        result = 1;
    end
elseif(exitflag==5) % CPLEX success with numerical issues
    disp(strcat(['Exit flag ',num2str(exitflag),' in ISEMPTY']));
    norm_csi = x(1);
    if norm_csi <= 1
        result = 0;
    else
        result = 1;
    end    
else % Not success (treated as somehow unfeasible and the point z not belonging to the line zonotope Z)
    result = 1;
end


% First version: 26-07-2020  


end




