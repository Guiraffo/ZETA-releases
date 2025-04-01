function deltahat_u = LTI_find_upperbounds(Utilde,Ntilde,Ztilde,N,deltahat_g,imodel)
%LTI_FIND_UPPERBOUNDS find the upper bounds for the deltahat variable for
%                     each pair of linear models and given input bounds
%                     (solves eq (34) as a MILP for each q)
%
%   SYNTAX: deltahatupper = LTI_FIND_UPPERBOUNDS(Utilde,Ntilde,Ztilde)
%
%   INPUTS
%      Utilde: augmented set of admissible inputs (Polyhedron)
%      Ntilde: cell array of matrices Ntilde
%      Ztilde: cell array of zonotopes Ztilde
%           N: separation horizon
%  deltahat_g: initial guess for the upper bound
%      imodel: model to be suppressed from isolation (empty if no model should be suppressed)
%              
%   OUTPUT
%  deltahat_u: upper bound for the deltahat variable
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

OPTIONS.Display = 'off';

nof_models = size(Ntilde,2);
nof_inputs = size(Utilde.H,2)/N;


% Guess step multiplier and tolerance value
increase_ratio = 10;
epsilon = 1e-5;

Gq = cell(nof_models,nof_models);
cq = cell(nof_models,nof_models);
deltahat_u = cell(nof_models,nof_models);
deltahat_g = repmat({deltahat_g},nof_models,nof_models);

if(isempty(imodel))
    imodel = 0;
end

for i=1:nof_models-1
    for j=i+1:nof_models
        
        if((i~=imodel)&&(j~=imodel))        

            Gq{i,j} = Ztilde{i,j}.G;
            cq{i,j} = Ztilde{i,j}.c;

            % MILP for a pair of models
            
            space_dimension = size(Gq{i,j},1);
            nof_generators = size(Gq{i,j},2);

            deltahat_u{i,j} = Inf;

            while(deltahat_u{i,j} + epsilon > deltahat_g{i,j})

                deltahat_g{i,j} = increase_ratio*deltahat_g{i,j};

                f = [zeros(nof_inputs*N,1); 1; zeros(1 + 5*nof_generators + space_dimension, 1)];

                Aineq = [  Utilde.H, zeros(size(Utilde.H,1), 1), zeros(size(Utilde.H,1), 1), zeros(size(Utilde.H,1), nof_generators), zeros(size(Utilde.H,1), space_dimension), zeros(size(Utilde.H,1), nof_generators), zeros(size(Utilde.H,1), nof_generators), zeros(size(Utilde.H,1), nof_generators), zeros(size(Utilde.H,1), nof_generators);
                         zeros(1, nof_inputs*N),  1, 0, zeros(1, nof_generators), zeros(1, space_dimension), zeros(1, nof_generators), zeros(1, nof_generators), zeros(1, nof_generators), zeros(1, nof_generators);            
                         zeros(1, nof_inputs*N), -1, 1, zeros(1, nof_generators), zeros(1, space_dimension), zeros(1, nof_generators), zeros(1, nof_generators), zeros(1, nof_generators), zeros(1, nof_generators);
                         zeros(nof_generators, nof_inputs*N),  zeros(nof_generators,1), -ones(nof_generators,1),      eye(nof_generators), zeros(nof_generators, space_dimension),  zeros(nof_generators, nof_generators), zeros(nof_generators, nof_generators), zeros(nof_generators, nof_generators), zeros(nof_generators, nof_generators);
                         zeros(nof_generators, nof_inputs*N),  zeros(nof_generators,1), -ones(nof_generators,1),     -eye(nof_generators), zeros(nof_generators, space_dimension),  zeros(nof_generators, nof_generators), zeros(nof_generators, nof_generators), zeros(nof_generators, nof_generators), zeros(nof_generators, nof_generators);
                         zeros(nof_generators, nof_inputs*N),  zeros(nof_generators,1), zeros(nof_generators,1),    zeros(nof_generators), zeros(nof_generators, space_dimension),  -eye(nof_generators), zeros(nof_generators), zeros(nof_generators), zeros(nof_generators);
                         zeros(nof_generators, nof_inputs*N),  zeros(nof_generators,1), zeros(nof_generators,1),    zeros(nof_generators), zeros(nof_generators, space_dimension), zeros(nof_generators),  -eye(nof_generators), zeros(nof_generators), zeros(nof_generators);
                         zeros(nof_generators, nof_inputs*N),  zeros(nof_generators,1), zeros(nof_generators,1),    zeros(nof_generators), zeros(nof_generators, space_dimension),   eye(nof_generators), zeros(nof_generators),  -eye(nof_generators), zeros(nof_generators);
                         zeros(nof_generators, nof_inputs*N),  zeros(nof_generators,1), zeros(nof_generators,1),    zeros(nof_generators), zeros(nof_generators, space_dimension), zeros(nof_generators),   eye(nof_generators), zeros(nof_generators),  -eye(nof_generators);                  
                         zeros(nof_generators, nof_inputs*N),   ones(nof_generators,1), zeros(nof_generators,1),     -eye(nof_generators), zeros(nof_generators, space_dimension), zeros(nof_generators), zeros(nof_generators), 2*(1+deltahat_g{i,j})*eye(nof_generators),                     zeros(nof_generators);
                         zeros(nof_generators, nof_inputs*N),   ones(nof_generators,1), zeros(nof_generators,1),      eye(nof_generators), zeros(nof_generators, space_dimension), zeros(nof_generators), zeros(nof_generators),                     zeros(nof_generators), 2*(1+deltahat_g{i,j})*eye(nof_generators);
                         zeros(nof_generators, nof_inputs*N),  -ones(nof_generators,1), zeros(nof_generators,1),      eye(nof_generators), zeros(nof_generators, space_dimension), zeros(nof_generators), zeros(nof_generators),  zeros(nof_generators), zeros(nof_generators);
                         zeros(nof_generators, nof_inputs*N),  -ones(nof_generators,1), zeros(nof_generators,1),     -eye(nof_generators), zeros(nof_generators, space_dimension), zeros(nof_generators), zeros(nof_generators),  zeros(nof_generators), zeros(nof_generators)];

                bineq = [               Utilde.k;
                                 deltahat_g{i,j};            
                                               1;
                         zeros(nof_generators,1);
                         zeros(nof_generators,1);
                         zeros(nof_generators,1);
                         zeros(nof_generators,1);
                         zeros(nof_generators,1);
                         zeros(nof_generators,1);      
                         (-1 + 2*(1+deltahat_g{i,j}))*ones(nof_generators,1);
                         (-1 + 2*(1+deltahat_g{i,j}))*ones(nof_generators,1);
                          ones(nof_generators,1);
                          ones(nof_generators,1)];

                Aeq = [     Ntilde{i,j}, zeros(space_dimension,1), zeros(space_dimension,1), -Gq{i,j}, zeros(space_dimension), zeros(space_dimension,nof_generators), zeros(space_dimension,nof_generators), zeros(space_dimension,nof_generators), zeros(space_dimension,nof_generators);
                       zeros(nof_generators, nof_inputs*N), zeros(nof_generators,1), zeros(nof_generators,1), zeros(nof_generators), Gq{i,j}.', -eye(nof_generators), eye(nof_generators), zeros(nof_generators), zeros(nof_generators);
                       zeros(1, nof_inputs*N), 0, 0, zeros(1, nof_generators), zeros(1, space_dimension), ones(1, nof_generators), ones(1, nof_generators), zeros(1,nof_generators), zeros(1,nof_generators)];

                beq = [                 cq{i,j};
                        zeros(nof_generators,1);
                                              1];

                ctype = strcat([repmat('C',1,nof_inputs*N+1+1+nof_generators+space_dimension+2*nof_generators), repmat('B',1,2*nof_generators)]);
                
                [x,~,exitflag] = optim.solvemilp(-f,Aineq,bineq,Aeq,beq,ctype,OPTIONS);

                if(ischar(exitflag)) % Gurobi
                    if(strcmp(exitflag,'OPTIMAL')) % Success
                        deltahat_u{i,j} = x(nof_inputs*N+1,1);
                    else % Not success 
                        error(strcat(['Exit flag ',exitflag,' in LTI_SEPARATING_INPUTS']));
                    end    
                elseif(exitflag==1) % CPLEX: Success
                    deltahat_u{i,j} = x(nof_inputs*N+1,1);  
                    %disp(strcat(['Exit flag ',num2str(exitflag),' in LTI_SEPARATING_INPUTS']));
                elseif(exitflag==5) % CPLEX: Solution with numerical issues
                    deltahat_u{i,j} = x(nof_inputs*N+1,1); 
                    disp(strcat(['Exit flag ',num2str(exitflag),' in LTI_SEPARATING_INPUTS']));
                else % Not success 
                    error(strcat(['Exit flag ',num2str(exitflag),' in LTI_SEPARATING_INPUTS']));
                end  

            end
        
        end

    end
end


%Q = nchoosek(nof_models,2);


% Revison 07-12-2018: Updated to U described by a Polyhedron instead of a box
% First version: 31-07-2018

end
