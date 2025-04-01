function Z = sec(X)
%SEC returns the secant of an interval
%
%   SYNTAX: Z = SEC(X)
%
%   INPUTS
%           Z: the input interval
%
%   OUTPUT
%           Z: the secant of the input interval

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

% Input check
xmin = -pi/2;
xmax = pi/2;
xL = inf(X);
xU = sup(X);

% Check if input is not contained in [-pi/2,pi/2]
if((xL<=xmin)||(xL>=xmax)||(xU<=xmin)||(xU>=xmax))
    error(' The secant function is currently implemented only in the interval [-pi/2,pi/2].');
end

if(xU<0) % Negative input
    Z = Interval(sec(xU),sec(xL));
elseif(xL>0) % Positive input
    Z = Interval(sec(xL),sec(xU));
else % If input contains zero
    Z = Interval(0,max(sec(xL),sec(xU)));
end

% 20-05-2024: First version
    
end