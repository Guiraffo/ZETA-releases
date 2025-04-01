function [x,fval,exitflag,output] = solvemilp(f,Aineq,bineq,Aeq,beq,ctype,OPTIONS)
%SOLVEMILP relays a mixed-integer linear program to a solver based on this
%          toolbox setup parameters
%
%   SYNTAX: [x,fval,exitflag,output] = SOLVEMILP(f,Aineq,bineq,Aeq,beq,LB,UB,OPTIONS)
%

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

switch toolsettings.solvermilp
    case 'intlinprog'
        persistent OPTIONS_;
        if(isempty(OPTIONS_))
            OPTIONS_ = optimoptions('intlinprog');
        end
        OPTIONS_.Display = OPTIONS.Display;
        [x,fval,exitflag,output] = intlinprog(f,find(ctype=='B'),Aineq,bineq,Aeq,beq,[],[],OPTIONS_);
    case 'gurobi'
        [x,fval,exitflag,output] = optim.gurobimilp(f,Aineq,bineq,Aeq,beq,ctype,OPTIONS.Display);
    otherwise
        error('Invalid MILP solver parameter.')
end
        

end
