function [Znew,P_vert,hz,lambda] = propagate(Z,f,OPTIONS)
%PROPAGATE performs set propagation based on DC programming principles.
%          (Alamo et al. 2005, de Paula et al., 2024). Internal use only. 

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

switch OPTIONS.vertices
    case 'hull'
        
        P = Zonotope(intervalhull(Z));
        P_vert = util.partope_vertices(P);
        hz = P.c;
             
    case 'partope'
        
        P = partopebound(Z);
        P_vert = util.partope_vertices(P);
        hz = P.c;
        
    otherwise
        error('Invalid algorithm option in DC programming propagate.');
end

switch OPTIONS.decomposition
    case 'exact'
        lambda = [];
    case 'aBB'
        lambda = dcprog.getlambda_alphaBB(f,intervalhull(Z));
    otherwise
        error('Invalid decomposition mode.');
end

[e_minus,e_plus] = dcprog.get_error_bounds(f,P_vert,hz,OPTIONS,lambda);
Rem = Zonotope(0.5*(e_plus+e_minus), 0.5*diag(e_plus-e_minus));

f_ath = f.eval(hz);
f_Jacob_ath = f.Jacob(hz);

PSI_linear = (f_ath - f_Jacob_ath*hz) + f_Jacob_ath*Z;
Znew = PSI_linear + Rem;

end

