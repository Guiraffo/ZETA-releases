function [Radius,Radius1norm,Hull] = radius(Z)
%RADIUS returns the radius of a zonotope (see below for different metrics)
%
%   SYNTAX: [Radius,Radius1norm,Hull] = RADIUS(Z)
%
%   INPUTS
%           Z: zonotope
%              
%   OUTPUT
%         Radius: the radius of Z (half the maximum width of the interval
%                 hull)
%    Radius1norm: the 1-norm version of the radius of Z
%           Hull: the interval hull of Z (for convenience use)

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

end


