function Znew = reduction(Z,ng_max,nc_max,ng_prior)
%REDUCTION performs complexity reduction of a constrained zonotope
%
%   SYNTAX: Znew = REDUCTION(Z,ng_max,nc_max)
%           Znew = REDUCTION(Z,ng_max,nc_max,ng_prior)
%
%   INPUTS
%           Z: constrained zonotope as a structure
%      ng_max: desired number of generators (empty if no generator reduction)
%      nc_max: desired number of constraints (empty if no constraint elimination)
%    ng_prior: desired number of generators for prior generator reduction
%              of CZs which have a huge amount of generators
%              (hundreds, thousands, leave it empty if not desired,
%              default is empty)
%
%   OUTPUT
%        Znew: the new constrained zonotope

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

% Performs constraint elimination then generator reduction


if(nargin<3)
    error('Not enough input arguments.');
elseif(nargin==3)
    ng_prior = [];
end


if(~isempty(ng_prior))
    Z = reduction_gen(Z,ng_prior);
end
if(~isempty(nc_max))
    Z = reduction_con(Z,nc_max);
end
if(~isempty(ng_max))
    Z = reduction_gen(Z,ng_max);
end
Znew = Z;

end

function Znew = reduction_con(Z,nc_max)
%REDUCTION_CON uses libCZon.ScaleDualize for constraint elimination as
%              described in the 2020 Automatica paper
%
%   SYNTAX: Znew = REDUCTION_CON(Z,nc_max)
%
%   INPUTS
%           Z: constrained zonotope as a structure
%      nc_max: desired number of constraints
%
%   OUTPUT
%        Znew: the new constrained zonotope


% Wrapper to iterated use of libCZon.ScaleDualize

nc_diff = size(Z.A,1) - nc_max;

if(nc_diff <= 0)
    Znew = Z;
    return;
end

Z_ = libCZon.Wrapto(Z);

for i=1:nc_diff
    
    Z_ = libCZon.ScaleDualize(Z_);
    
    % Stop condition due to the possibility of removing null constraints inside libCZon_ScaleDualize3
    if(size(Z_{4},1)<=nc_max)
        break;
    end
    
end

Znew = libCZon.Wrapfrom(Z_);

end

function Znew = reduction_gen(Z,ng_max)
%REDUCTION_GEN performs generator reduction of a constrained zonotope
%
%   SYNTAX: Znew = REDUCTION_GEN(Z,ng_max)
%
%   INPUTS
%           Z: constrained zonotope as a structure
%      ng_max: desired number of generators
%
%   OUTPUT
%        Znew: the new constrained zonotope


% Wrapper to libZon.ReduceOrderChisci for generator reduction of CZ


ng_diff = size(Z.G,2) - ng_max;

if(ng_diff <= 0)
    Znew = Z;
    return;
end


zo = ng_max/(size(Z.G,1) + size(Z.A,1));

[Zlift,dimension] = lift(Z);

Zlift_    = libZon.Wrapto(Zlift);
Zliftred_ = libZon.ReduceOrderChisci(Zlift_, zo);
Zliftred  = libZon.Wrapfrom(Zliftred_);

Znew = unlift(Zliftred,dimension);

end


