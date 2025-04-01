function [h,Z,Zbar] = chooseh(Z,Zhull,Jacob,heuristic)
%CHOOSEH returns the approximation point h in the interval hull of Z for a
%        given interval matrix and heuristic (for CZMV and CZFO methods)
%
%   SYNTAX: h = CHOOSEH(Z,Zhull,Jacob,heuristic)
%
%   INPUTS
%           Z: the contrained zonotope
%       Zhull: the interval hull of Z (interval object)
%       Jacob: interval matrix enclosing the range of the respective
%              Jacobian (required by 'mindiamm' only)
%   heuristic:    'midhull': returns the center of the interval hull
%              'distcenter': returns the closest point to c in the interval
%                            hull
%                'mindiamm': returns the point in interval hull that
%                            minimizes the diameter of m (CZ-inclusion)
%
%   OUTPUT
%           h: the real vector h in the interval hull of Z according to the
%              chosen heuristic
%           Z: the associated equivalent CG-rep of Z 

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


if(strcmp(heuristic,'midhull'))

    h = mid(Zhull);
    
elseif(strcmp(heuristic,'distcenter'))
    
        
    if(isinside(Zhull,Z.c)) % If c belongs to the interval hull, use c

        h = Z.c;

    else % If not, use the closest point to c, according to 1-norm distance
        
        h = closest(Z,Z.c,'hull',Zhull);

    end

elseif(strcmp(heuristic,'mindiamm'))

    nof_constraints = size(Z.A,1);
    Gbar = cell(nof_constraints,1);
    Abar = cell(nof_constraints,1);
    bbar = cell(nof_constraints,1);

    Gtilde = cell(nof_constraints,1);
    Atilde = cell(nof_constraints,1);
    btilde = cell(nof_constraints,1);        

    xi_m = cell(nof_constraints,1);
    xi_r = cell(nof_constraints,1);
    Lambda_G = cell(nof_constraints,1);

    Zred_ = libCZon.Wrapto(Z);
    j = 1;
    while(j<=nof_constraints)
        [Zred_,~,~,flag,xi_elim,xi_m{j},xi_r{j},Gbar{j},Abar{j},bbar{j},pivot_row] = libCZon.ScaleDualize(Zred_);

        % Take into account constraints that were removed because were null (indicated by flag = -9)
        if(flag==-9) % If null constraints were removed after rescaling
            nof_removednullconstr = length(pivot_row)-1;

            Gtilde{j} = Gbar{j}*diag(xi_r{j});
            Atilde{j} = Abar{j}*diag(xi_r{j});
            btilde{j} = bbar{j} - Abar{j}*xi_m{j};  

            Lambda_G{j} = zeros(size(Gbar{j},1),size(Abar{j},1)); % No generators were removed for this one

            for m=1:nof_removednullconstr-1 % Gauss elimination and rescaling were performed only before removing the first null constraint

                Gbar{j+m} = Gbar{j};
                Abar{j+m} = Abar{j};
                bbar{j+m} = bbar{j};

                Gtilde{j+m} = Gtilde{j};
                Atilde{j+m} = Atilde{j};
                btilde{j+m} = btilde{j};

                xi_r{j+m} = ones(size(xi_r{j},1),1);
                xi_m{j+m} = zeros(size(xi_m{j},1),1);

                Lambda_G{j+m} = zeros(size(Gbar{j},1),size(Abar{j},1)); % No generators were removed for these ones

            end

            if(pivot_row(end)~=-1) % If a non-null constraint was removed besides null constraints. -1 indicates that only null constraints were removed

                % Gauss elimination and rescaling were not performed

                Gbar{j+nof_removednullconstr} = Gbar{j};
                Abar{j+nof_removednullconstr} = Abar{j};
                bbar{j+nof_removednullconstr} = bbar{j};

                Gtilde{j+nof_removednullconstr} = Gtilde{j};
                Atilde{j+nof_removednullconstr} = Atilde{j};
                btilde{j+nof_removednullconstr} = btilde{j};

                xi_r{j+nof_removednullconstr} = ones(size(xi_r{j},1),1);
                xi_m{j+nof_removednullconstr} = zeros(size(xi_m{j},1),1);


                Ej1 = zeros(size(Gbar{j},2),size(Abar{j},1));
                Ej1(xi_elim,pivot_row(end)) = 1;
                Lambda_G{j+nof_removednullconstr} = Gtilde{j+nof_removednullconstr}*Ej1/Atilde{j+nof_removednullconstr}(pivot_row(end),xi_elim); % Corrected pivot row                

                j = j + nof_removednullconstr + 1;  

            else

                j = j + nof_removednullconstr;

            end

        else

            Gtilde{j} = Gbar{j}*diag(xi_r{j});
            Atilde{j} = Abar{j}*diag(xi_r{j});
            btilde{j} = bbar{j} - Abar{j}*xi_m{j};

            Ej1 = zeros(size(Gbar{j},2),size(Abar{j},1));
            Ej1(xi_elim,pivot_row) = 1;
            Lambda_G{j} = Gtilde{j}*Ej1/Atilde{j}(pivot_row,xi_elim); % Corrected pivot row

            j = j + 1;

        end

    end

    hstar = Z.c;
    for j=1:nof_constraints
        hstar = hstar + Gbar{j}*xi_m{j} + Lambda_G{j}*btilde{j};
    end        

    if(isinside(Zhull,hstar)) % If hstar belongs to the interval hull, use hstar

        h = hstar;

    else % If not, use the point h in the interval hull that minimizes ||diam(m)||_1

        space_dimension = size(Z.G,1);

        Theta = diag(sum(diam(Jacob),1));

        f = [ones(space_dimension,1); zeros(space_dimension,1)]; % Ordering: l, gamma

        A = [-eye(space_dimension),  Theta; % 1-norm constraints 
             -eye(space_dimension), -Theta];

        b = [  Theta*hstar;
              -Theta*hstar];

        Aeq = [];
        beq = [];

        LB = [ zeros(space_dimension,1); inf(Zhull)];
        UB = [   Inf(space_dimension,1); sup(Zhull)];


        OPTIONS.Display = 'off';
        [x,~,exitflag] = optim.solvelp(f,A,b,Aeq,beq,LB,UB,OPTIONS);

        if(ischar(exitflag)) % Gurobi
            if(~strcmp(exitflag,'OPTIMAL')) % Success    norm_csi = x(1);        
                error(strcat(['Error in CZonotope/chooseh "mindiamm": LP exitflag is ',num2str(exitflag)]));
            end
        elseif(exitflag~=1) % Linprog
            error(strcat(['Error in CZonotope/chooseh "mindiamm": LP exitflag is ',num2str(exitflag)]));
        end
 
        h = x(space_dimension+1:end);   
           
    end
    
    % Return Zbar as well for use in the CZ-inclusion
    Zbar = libCZon.Wrapfrom(Zred_);
    Zbar.c = Zbar.c - h; % Displace the center by -h
    
else
    error('Invalid heuristic for h in CZonotope/chooseh.');
end    



end
