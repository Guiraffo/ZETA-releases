function val = faultdiagopts(functional,R,N,max_order,epsilon)
%FAULTDIAGOPTS returns an OPTIONS structure for fault diagnosis methods
%
%   INPUTS
%   functional: use a 'quadratic' or a 'linear' cost function
%            R: cost function weighting matrix
%            N: separation horizon
%    max_order: max allowed order for the zonotopes
%      epsilon: minimum separation threshold

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

val.functional = functional;
val.R = R;
val.N = N;
val.max_order = max_order;
val.epsilon = epsilon;
val.deltahat_guess = 100; % Default value

end

