function val = polyrelax_approxmode(newmode)
% Set/get approximation mode for the linearized convex/concave enclosure
% Modes:
% 0 = 3 points (zL,zmid,zU);
%   - Default mode. Higher H-rep complexity. THIS IS THE ONLY FULLY
%     SUPPORTED MODE FOR NOW.
% 1 = 2 points (zL,zU)
%   - Reduced H-rep complexity. Potentially poor enclosure.
% 2 = one point at 1/3, another at 2/3 ((1/3)*(zL+zU), (2/3)*(zL+zU))
%   - Same complexity as Mode 1. However, this mode requires the
%     intersection with interval Z at the end of the algorithm (as it is
%     already done in the CZ reachability/estimation algorithm) to properly
%     clip the enclosure, to give a better result than Mode 1.

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

mlock; % So clear doesnt reset the static variable
persistent approxmode_;
if(isempty(approxmode_))
    approxmode_ = 0; % Default mode
end
if nargin
    approxmode_ = newmode;
end
val = approxmode_;  
    
% 21-05-2024: first version    
    
end  

