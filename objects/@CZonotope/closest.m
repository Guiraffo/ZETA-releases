function zstar = closest(Z,h,bound,Zhull)
%CLOSEST returns the closest point inside a given constrained zonotope
%        with respect to a desired point h (distance in the 1-norm sense)
%
%   SYNTAX: zstar = CLOSEST(Z,h,bound,Zhull)
%
%   INPUTS
%           Z: constrained zonotope object
%           h: point of interest
%       bound: 'set': find the closest point in Z
%             'hull': find the closest point in the interval hull of Z
%       Zhull: interval hull of Z (leave it empty if not required)
%
%   OUTPUT
%       zstar: the closest point to h

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


if(strcmp(bound,'set'))



    % OPTIONS
    OPTIONS.Display = 'off';

    space_dimension = size(Z.G,1);
    nof_generators = size(Z.G,2);
    nof_constraints = size(Z.A,1);

    f = [ones(space_dimension,1); 0; zeros(nof_generators,1)]; % Ordering: l, t, csi

    A = [-eye(space_dimension), zeros(space_dimension,1), -Z.G; % 1-norm constraints 
         -eye(space_dimension), zeros(space_dimension,1),  Z.G;
         -eye(space_dimension), zeros(space_dimension,1), zeros(space_dimension,nof_generators);
         zeros(nof_generators,space_dimension), -ones(nof_generators,1),  eye(nof_generators); % Inf-norm constraints
         zeros(nof_generators,space_dimension), -ones(nof_generators,1), -eye(nof_generators);
                      zeros(1,space_dimension),                       1,  zeros(1,nof_generators); 
                      zeros(1,space_dimension),                      -1,  zeros(1,nof_generators)];

    b = [Z.c-h;-(Z.c-h);zeros(space_dimension + 2*nof_generators,1); 1; 0];

    Aeq = [zeros(nof_constraints,space_dimension), zeros(nof_constraints,1), Z.A];
    beq = Z.b;

    [x,~,exitflag] = optim.solvelp(f,A,b,Aeq,beq,[],[],OPTIONS);

    zstar = Z.c + Z.G*x(space_dimension+2:end); 
    
    
elseif(strcmp(bound,'hull'))
    
    if(nargin==3)
        Zhull = [];
    end
    if(isempty(Zhull))
        Zhull = intervalhull(Z);
    end
    
    % OPTIONS
    OPTIONS.Display = 'off';

    space_dimension = size(Z.G,1);

    f = [ones(space_dimension,1); zeros(space_dimension,1)]; % Ordering: l, z

    A = [-eye(space_dimension), -eye(space_dimension); % 1-norm constraints 
         -eye(space_dimension),  eye(space_dimension)];

    b = [-h
          h];

    Aeq = [];
    beq = [];
    
    LB = [-Inf(space_dimension,1); inf(Zhull)];
    UB = [ Inf(space_dimension,1); sup(Zhull)];

    [x,~,exitflag] = optim.solvelp(f,A,b,Aeq,beq,LB,UB,OPTIONS);
    
    zstar = x(space_dimension+1:end);
    
else
    
    error('Invalid "bound" argument in CZonotope/closest.'); 
    
end
    
end




