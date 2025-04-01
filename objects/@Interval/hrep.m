function [H,k] = hrep(B)
%HREP converts an interval into a polytope in halfspace representation
%
%   SYNTAX: V = HREP(B)
%
%   INPUTS
%           B: box as an interval column vector
%
%   OUTPUT
%         H,k: interval in halfspace representation

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

dimension = size(B,1);
if(size(B,2)>1)
    error('Invalid input dimensions in Interval/hrep.');
end

H = [ eye(dimension);
     -eye(dimension)];
k = [ sup(B);
     -inf(B)];

end