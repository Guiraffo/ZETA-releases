function [Radius,Radius1norm,Hull,HullVol] = radius(Z)
%RADIUS returns the radius of a constrained zonotope according to different
%       metrics
%
%   SYNTAX: [Radius] = RADIUS(Z)
%           [Radius,Radius1norm] = RADIUS(Z)
%           [Radius,Radius1norm,Hull] = RADIUS(Z)
%           [Radius,Radius1norm,Hull,HullVol] = RADIUS(Z)
%
%   INPUTS
%           Z: constrained zonotope as a structure object
%              
%   OUTPUT
%         Radius: half the width of the interval hull
%    Radius1norm: the 1-norm version of the radius
%           Hull: the interval hull of Z (for convenience use)
%        HullVol: the volume of the interval hull of Z (for convenience use)

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


Hull = intervalhull(Z);
Radius = max(rad(Hull));
Radius1norm = sum(rad(Hull));
HullVol = prod(diam(Hull));


end


