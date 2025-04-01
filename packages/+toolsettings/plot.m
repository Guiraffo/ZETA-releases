function val = plot(setplot)
% Set/get plot method
% Allowed values for setplot:
% 'yalmip'  = use yalmip
% 'mpt-H'   = use MPT halfspace representation, if available
% 'mpt-V'   = use MPT vertex representation (will use mtp-H if unavailable)
% Default is 'yalmip'

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
persistent plotmethod_;
if(isempty(plotmethod_))
    plotmethod_ = 'yalmip'; % Default plot method
end
if nargin
    plotmethod_ = setplot;
end
val = plotmethod_;  
        
% 03-10-2024: first version

end

