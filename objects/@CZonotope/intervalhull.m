function Zeta = intervalhull(Z)
%INTERVALHULL returns the smallest box containing a constrained zonotope
%             using linear programming
%
%   SYNTAX: Zeta = INTERVALHULL(Z)
%
%   INPUTS
%           Z: constrained zonotope object
%              
%   OUTPUT
%        Zeta: the smallest interval containing the constrained zonotope 

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
   
space_dimension = size(Z.G,1);
nof_generators = size(Z.G,2);
nof_constraints = size(Z.A,1);

zetalower = zeros(space_dimension,1);
zetaupper = zeros(space_dimension,1);

for j=1:space_dimension

    f = [0; Z.G(j,:).'];
    A = [                      1,  zeros(1,nof_generators);
                              -1,  zeros(1,nof_generators);
         -ones(nof_generators,1),      eye(nof_generators);
         -ones(nof_generators,1),     -eye(nof_generators)];
    b = [1; 0; zeros(2*nof_generators,1)];
    Aeq = [zeros(nof_constraints,1), Z.A];
    beq = Z.b;

    LB = [0; -ones(nof_generators,1)];
    UB = [1;  ones(nof_generators,1)];

    
    [~,zmin,exitflag] = optim.solvelp(f,A,b,Aeq,beq,LB,UB,OPTIONS);  zetalower(j) =  zmin; 

    if(ischar(exitflag)) % Gurobi
        if(~strcmp(exitflag,'OPTIMAL')) % Success
            error(strcat(['Error in CZonotope_intervalhull. Exitflag = ',exitflag]));
        end          
    elseif(exitflag~=1) % Linprog
        error(strcat(['Error in CZonotope/intervalhull. Exitflag = ',num2str(exitflag)]));
    end     

    [~,zmax,exitflag] = optim.solvelp(-f,A,b,Aeq,beq,LB,UB,OPTIONS); zetaupper(j) = -zmax;        

    if(ischar(exitflag)) % Gurobi
        if(~strcmp(exitflag,'OPTIMAL')) % Success
            error(strcat(['Error in CZonotope/intervalhull. Exitflag = ',exitflag]));
        end          
    elseif(exitflag~=1) % Linprog
        error(strcat(['Error in CZonotope/intervalhull. Exitflag = ',num2str(exitflag)]));
    end


end

zetalower = zetalower + Z.c;
zetaupper = zetaupper + Z.c;

switch toolsettings.IA_class
    case 'intval'
        Zeta = infsup(zetalower,zetaupper);
    case 'Interval'
        Zeta = Interval(zetalower,zetaupper);
    otherwise
        error('Invalid IA class in CZonotope/intervalhull');
end   
    
end



