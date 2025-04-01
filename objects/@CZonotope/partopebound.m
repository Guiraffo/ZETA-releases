function P = partopebound(Z,alg)
%PARTOPEBOUND returns a parallelotope bounding a constrained zonotope
%
%   SYNTAX: P = PARTOPEBOUND(Z)
%
%   INPUTS
%           Z: constrained zonotope object
%              
%   OUTPUT
%           P: parallelotope bounding the constrained zonotope (in CG-rep)

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

%dim = size(Z.G,1);
%nof_g = size(Z.G,2);
%nof_constraints = size(Z.A,1);


if(nargin==1)
    alg = 'rescale';
end

switch alg
        
    case 'red'
        
        % Using complexity reduction: eliminate all constraints, then
        % reduce the number of generators to a paralelotope
        
        P = reduction(Z,Z.dim,0);
        
    case 'rescale'
        
        % Rescale-based method, Proposition 3 in de Paula et al (2024)
        
        
        P_outer = partopebound(Z,'red');
        
        Gc = P_outer.G;
        cc = P_outer.c;
        Gz = Z.G;
        cz = Z.c;
        Az = Z.A;
        bz = Z.b;
        
        ng_c = size(Gc,2);        
        ng_z = size(Gz,2);
        nc_z = size(Az,1);
        
        % Options
        OPTIONS.Display = 'off';
        
        % Decision variables
        % Generators of P_outer and Z, xic and xiz        
        
        % Constraints
        % Gc*xic + cc = Gz*xiz + cz
        % xic \in [-1,1], xiz \in [-1,1]
        
        % Inf norm constraints
        Aineq = [       -eye(ng_c), zeros(ng_c,ng_z);
                         eye(ng_c), zeros(ng_c,ng_z);
                  zeros(ng_z,ng_c),       -eye(ng_z);
                  zeros(ng_z,ng_c),        eye(ng_z)];    
              
        bineq = ones(2*(ng_z+ng_c),1);
        
        % Equality constraints
        
        Aeq = [              Gc, -Gz;
               zeros(nc_z,ng_c),  Az];
        beq = [cz - cc;
                    bz];
                
        zetaL = zeros(ng_c,1);
        zetaU = zeros(ng_c,1);
                
        for j=1:ng_c
            
            fobj = zeros(ng_c+ng_z,1);
            fobj(j) = 1;
            
            [~,foptimal,exitflag,~] = optim.solvelp( fobj,Aineq,bineq,Aeq,beq,[],[],OPTIONS); zetaL(j) =  foptimal;
            [~,foptimal,exitflag,~] = optim.solvelp(-fobj,Aineq,bineq,Aeq,beq,[],[],OPTIONS); zetaU(j) = -foptimal;
            
        end
            
        P = CZonotope(cc + Gc*(0.5*(zetaU+zetaL)), Gc*diag(0.5*(zetaU-zetaL)));
        
    otherwise
        
        error('Invalid algorithm in CZonotope/partopebound');
        
end

end



