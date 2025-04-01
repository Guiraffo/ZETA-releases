function Zeta = intervalhull(Z)
%INTERVALHULL returns the smallest box containing a line zonotope using
%             linear programming. Can be infinite if the set is unbounded.
%
%   SYNTAX: Zeta = INTERVALHULL(Z)
%
%   INPUTS
%           Z: line zonotope object
%              
%   OUTPUT
%        Zeta: the smallest interval containing the line zonotope. Limits
%              will be [-Inf, Inf] if the line zonotope is unbounded.

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


persistent OPTIONS; % Helps computational time
if(isempty(OPTIONS))
    OPTIONS.Display = 'off';
end
    
space_dim = size(Z.G,1);
nof_gens = size(Z.G,2);
nof_cons = size(Z.A,1);
nof_lins = size(Z.M,2);

zetalower = zeros(space_dim,1);
zetaupper = zeros(space_dim,1);

for j=1:space_dim

    f = [0; Z.G(j,:).'; Z.M(j,:).'];
    A = [                  1,  zeros(1,nof_gens), zeros(1,nof_lins);
                          -1,  zeros(1,nof_gens), zeros(1,nof_lins);
           -ones(nof_gens,1),      eye(nof_gens), zeros(nof_gens,nof_lins);
           -ones(nof_gens,1),     -eye(nof_gens), zeros(nof_gens,nof_lins)];
    b = [1; 0; zeros(2*nof_gens,1)];
    Aeq = [zeros(nof_cons,1), Z.A, Z.S];
    beq = Z.b;

    LB = [0; -ones(nof_gens,1); -Inf(nof_lins,1)];
    UB = [1;  ones(nof_gens,1);  Inf(nof_lins,1)];


    % Lower bound
    [~,zmin,exitflag] = optim.solvelp(f,A,b,Aeq,beq,LB,UB,OPTIONS);
    if(ischar(exitflag)) % Gurobi
        if(strcmp(exitflag,'OPTIMAL')) % Gurobi success
            zetalower(j) =  zmin;
        elseif(strcmp(exitflag,'INF_OR_UNBD')) % Gurobi unbounded
            zetalower(j) = -Inf;
        else % Gurobi failed
            error(strcat(['Error in LZonotope/intervalhull. Exitflag = ',exitflag]));
        end      
    elseif(exitflag==1) % Linprog success
        zetalower(j) = zmin;
    elseif(exitflag==-3) % Linprog unbounded
        zetalower(j) = -Inf;
    else % Linprog failed
        error(strcat(['Error in LZonotope/intervalhull. Exitflag = ',num2str(exitflag)]));
    end

    % Upper bound
    [~,zmax,exitflag] = optim.solvelp(-f,A,b,Aeq,beq,LB,UB,OPTIONS);    
    if(ischar(exitflag)) % Gurobi
        if(strcmp(exitflag,'OPTIMAL')) % Gurobi success
            zetaupper(j) = -zmax;
        elseif(strcmp(exitflag,'INF_OR_UNBD')) % Gurobi unbounded
            zetaupper(j) = Inf;
        else % Gurobi failed
            error(strcat(['Error in LZonotope/intervalhull. Exitflag = ',exitflag]));
        end      
    elseif(exitflag==1) % Linprog success
        zetaupper(j) = -zmax;
    elseif(exitflag==-3) % Linprog unbounded
        zetaupper(j) = Inf;
    else % Linprog failed
        error(strcat(['Error in LZonotope/intervalhull. Exitflag = ',num2str(exitflag)]));
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



