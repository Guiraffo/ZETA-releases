function Znew = reduction(Z,ng_max,method,W)
%REDUCTION performs generator reduction of a zonotope
%
%   SYNTAX: Znew = REDUCTION(Z,ng_max)
%           Znew = REDUCTION(Z,ng_max,method)
%           Znew = REDUCTION(Z,ng_max,method,W)
%
%   INPUTS
%           Z: zonotope as a structure object
%      ng_max: desired number of generators
%      method: - 'Combastel': uses Method 1 in Yang and Scott (2018)
%              - 'Girard': uses Method 2 in Yang and Scott (2018) 
%              - 'Chisci': uses Method 4 in Yang and Scott (2018)
%              - 'CombastelW': uses the weighted version of 'Combastel',
%                              presented in Combastel (2015)
%              Default method is 'Chisci'.
%           W: weighting matrix used in method 'Weighted'.
%
%   OUTPUT
%        Znew: the new zonotope

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


% Default method is 'Chisci'
if(nargin<2)
    error('Not enough inputs in Zonotope/reduction.')
elseif(nargin==2)
    method = 'Chisci';
else
    if(isempty(method))
        method = 'Chisci';
    end
end


ng_diff = size(Z.G,2) - ng_max;

if(ng_diff <= 0)
    Znew = Z;
    return;
end


switch method
    case 'Combastel'
        
        Znew_ = zonotope_reduction_combastel(Z,ng_max);        
        
    case 'Girard'
        
        zo = ng_max/(size(Z.G,1));
        Znew_ = libZon.Wrapfrom(libZon.ReduceOrderSimple(libZon.Wrapto(Z),zo));
        
    case 'Chisci'
        
        zo = ng_max/(size(Z.G,1));
        Znew_ = libZon.Wrapfrom(libZon.ReduceOrderChisci(libZon.Wrapto(Z),zo));
        
    case 'CombastelW'
        
        Znew_ = zonotope_reduction_combastel_weighted(Z,ng_max,W);        
        
    otherwise
        
        error('Invalid method in Zonotope/reduction.');
        
end

Znew = Zonotope(Znew_.c,Znew_.G);

% Revision 31-01-2020: added multiple methods (it was only 'Chisci')

end

function Z = zonotope_reduction_combastel(Zin,ngmax)
% Computes a lower-order zonotope that contains a given zonotope (Property
% 1 from Alamo et al, 2005)

if(size(Zin.G,2) <= ngmax)
    Z = Zin;
    return;
end

ng = size(Zin.G,2);
n = size(Zin.G,1);

Z.c = Zin.c;
Z.G = zeros(n,ngmax);

normvector = zeros(1,ng);
for i=1:ng    
    normvector(1,i) = norm(Zin.G(:,i),2); % 2-norm
end

[~,indexes] = sort(normvector,'descend');

Ghat = Zin.G(:,indexes(:)); % Sorted generators
GhatT = Ghat(:,1:ngmax-n); % Store bigger generators (according to 2-norm)

sumAbsGhat = sum(abs(Ghat(:,(ngmax-n+1):ng)),2); % Sum smaller generators
Q = diag(sumAbsGhat);

Z.G = [GhatT, Q];


end

function Z = zonotope_reduction_combastel_weighted(Zin,ngmax,W)
% Weighted version of the zonotope order reduction algorithm (Combastel, 
% 2015)

if(size(Zin.G,2) <= ngmax)
    Z = Zin;
    return;
end

ng = size(Zin.G,2);
n = size(Zin.G,1);

Z.c = Zin.c;
Z.G = zeros(n,ngmax);

normvector = zeros(1,ng);
for i=1:ng
    normvector(1,i) = sqrt(Zin.G(:,i).'*W*Zin.G(:,i)); % Weighted 2-norm
end

[~,indexes] = sort(normvector,'descend');

Ghat = Zin.G(:,indexes(:)); % Sorted generators
GhatT = Ghat(:,1:ngmax-n); % Store bigger generators (according to weighted 2-norm)

sumAbsGhat = sum(abs(Ghat(:,(ngmax-n+1):ng)),2); % Sum smaller generators
Q = diag(sumAbsGhat);

Z.G = [GhatT, Q];

end


