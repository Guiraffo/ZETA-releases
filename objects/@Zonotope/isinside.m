function result = isinside(Z,z)
%ISINSIDE verifies if a point belongs to a zonotope
%
%   SYNTAX: result = ISINSIDE(Z,z)
%
%   INPUTS
%           Z: zonotope object
%           z: point to be verified
%
%   OUTPUT
%      result: 0 if z does not belong to Z
%              1 if z belongs to Z

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

% Linear programming

space_dimension = Z.dim;
nof_generators = Z.ng;

% Linprog options
% OPTIONS = optimset('Display','off');

% OPTIONS
persistent OPTIONS; % Helps computational time
if(isempty(OPTIONS))
    OPTIONS.Display = 'off';
end

f = [1; zeros(nof_generators,1)];
A = [-ones(nof_generators,1),  eye(nof_generators);
     -ones(nof_generators,1), -eye(nof_generators)];
b = zeros(2*nof_generators,1);
Aeq = [zeros(space_dimension,1), Z.G];       
beq = z - Z.c;

% Solves the linear program
[x,~,exitflag] = optim.solvelp(f,A,b,Aeq,beq,[],[],OPTIONS);

if(ischar(exitflag)) % Gurobi
    if(strcmp(exitflag,'OPTIMAL')) % Success    norm_csi = x(1);
        norm_csi = x(1);
        if norm_csi <= 1%+eps
            result = 1;
        else
            result = 0;
        end
    else % Not success (treated as somehow unfeasible and the point z not belonging to the constrained zonotope Z)
        %disp(strcat(['Exit flag ',num2str(exitflag),' in CZONOTOPE_BELONGSTO']));
        result = 0; 
    end
elseif(exitflag==1) % CPLEX success
    norm_csi = x(1);
    if norm_csi <= 1%+eps
        result = 1;
    else
        result = 0;
    end
elseif(exitflag==5) % CPLEX success with numerical issues
    disp(strcat(['Exit flag ',num2str(exitflag),' in ZONOTOPE/ISINSIDE']));
    norm_csi = x(1);
    if norm_csi <= 1%+eps
        result = 1;
    else
        result = 0;
    end    
else % Not success (treated as somehow unfeasible and the point z not belonging to the constrained zonotope Z)
    %disp(strcat(['Exit flag ',num2str(exitflag),' in CZONOTOPE_BELONGSTO']));
    result = 0;
end

% Revision 13-12-2022: updated error treatment
% Revision 16-06-2020: renamed to szonotope_isinside
% 29-05-2018: changed solver to cplex

end




