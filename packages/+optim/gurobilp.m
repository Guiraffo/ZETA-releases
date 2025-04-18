function [x,fval,exitflag,output] = gurobilp(f,Aineq,bineq,Aeq,beq,LB,UB,Display)
%GUROBILP calls Gurobi to solve a linear program
%
%   SYNTAX: [x,fval,exitflag,output] = GUROBILP(f,Aineq,bineq,Aeq,beq,LB,UB,Display)
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

nof_variables = size(f,1);
nof_ineq = size(Aineq,1);
nof_eq = size(Aeq,1);

if(isempty(LB))
    LB = repmat(-Inf,nof_variables,1); % Default values are zero, need to overload
end
if(isempty(UB))
    UB = repmat( Inf,nof_variables,1);
end

model.obj = f.';
model.A = sparse([Aineq;Aeq]);
model.rhs = [bineq;beq];
model.sense = strcat([repmat('<',1,nof_ineq), repmat('=',1,nof_eq)]);
%model.vtype = ctype;
model.modelsense = 'min';
model.lb = LB;
model.ub = UB;

if(strcmp(Display,'on'))
    params.OutputFlag = 1;    
elseif(strcmp(Display,'off'))
    params.OutputFlag = 0;
else
    error('Invalid Display parameter in gurobilp');
end

%params.Method=0; % Uses primal simplex
%params.Method=1; % Uses dual simplex
%params.Method=2; % Uses barrier

results = gurobi(model,params);

if(strcmp(results.status,'OPTIMAL'))
    x = results.x;
    fval = results.objval;
    exitflag = results.status;
    output = results;
else
    x = [];
    fval = [];
    exitflag = results.status;
    output = [];    
end


end

