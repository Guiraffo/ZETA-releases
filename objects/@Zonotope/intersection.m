function Z = intersection(Zin,S,method)
%INTERSECTION computes a zonotope that contains the intersection
%                       of a zonotope and a strip
%
%   SYNTAX: Z = INTERSECTION(Zin,S)
%           Z = INTERSECTION(Zin,S,method)
%
%   INPUTS
%         Zin: zonotope object
%           S: strip object
%      method: 'seg': minimum segments (Alamo et al, 2005)
%            'Bravo': approximated minimum volume (Bravo et al, 2006)
%              (default is 'Bravo')
%
%   OUTPUT
%           Z: zonotope containing the intersection

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

if(nargin<2)
    error('Not enough input arguments.');
elseif(nargin==2)
    method = 'Bravo';
end
    
Zin = Zonotope(Zin);

if(~isa(S,'Strip'))
    error('Invalid inputs in Zonotope/intersection.')
end

if(strcmp(method, 'seg'))
    
    % Minimize the segments of the intersection zonotope

    Z = minimize_segments(Zin,S);
    
elseif(strcmp(method, 'Bravo'))
        
    % Minimize the volume of the intersection zonotope
        
    Z = minimize_volume(Zin,S);
        
else    
    
    error('Selected algorithm for intersection must be "seg" or "Bravo".');
        
end

end

% ----------- Minimization of the segments of the zonotope --------------
% Based on section 6.1 from 'Guaranteed State Estimation' by Alamo et al,
% 2005.
function Z = minimize_segments(Zin,S)

lambda = (Zin.G*(Zin.G.')*S.p)/((S.p.')*Zin.G*(Zin.G.')*S.p + S.sigma^2);

% Z.c = Zin.c + lambda*(S.d - S.p.'*Zin.c);
% Z.G = [(eye(size(Zin.G,1)) - lambda*(S.p.'))*Zin.G, S.sigma*lambda];

c = Zin.c + lambda*(S.d - S.p.'*Zin.c);
G = [(eye(size(Zin.G,1)) - lambda*(S.p.'))*Zin.G, S.sigma*lambda];

Z = Zonotope(c,G);

end


% --- Approximate minimum volume of the j-parametrized intersection ---
% Based on section III.C from 'Bounded Error Identification of Systems
% With Time-Varying Parameters' by Bravo et al, 2006.
function Z = minimize_volume(Zin,Strip)


    % Strip tightening
    SupportStrip = support_strip(Zin,Strip.p);
    Strip = intersection(Strip,SupportStrip);
  
    % Based on the original code in libCZon

    % Compute the intersection between the current zonotope and the tightened strip
    Qc = Strip.p.'*Zin.c;
    Qg = Strip.p.'*Zin.G;
    
    nof_generators = size(Zin.G,2);
    
    ZonCand = cell(nof_generators+1,1);
    Vol = zeros(nof_generators+1,1);
    ZonCand{1} = Zin;
    Vol(1) = abs(det(ZonCand{1}.G*ZonCand{1}.G')); % The first volume is the volume of the input zonotope
    
    % Compute ng candidate zonotopes
    for j=1:nof_generators
        if (abs(Qg(1,j))>1e-10)

            ZonCand{j+1} = Zonotope(Zin.c +( (Strip.d-Qc) / Qg(j) )*Zin.G(:,j),...
                                     Zin.G - kron( (Qg/Qg(j)) , Zin.G(:,j) ));
                               
            ZonCand{j+1}.G(:,j) = (Strip.sigma/Qg(j))*Zin.G(:,j);

        else
            
            ZonCand{j+1} = Zin;
            
        end
        
        Vol(j+1) = abs(det(ZonCand{j+1}.G*ZonCand{j+1}.G'));

    end

    %Determine the minimum volume zonotope
    [~, index] = sort(Vol);
    
    Z = ZonCand{index(1)}; % Gets the zonotope with smallest volume

end