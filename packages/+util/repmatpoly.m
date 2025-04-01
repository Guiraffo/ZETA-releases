function [Hrm,krm] = repmatpoly(H,k,N)
%REPMATPOLY returns the N-th iterative Cartesian product of a polyhedron
%           in halfspace representation {Hx<=k, Ax=b}
%
%   SYNTAX: INTERSECTPOLY(H,k,N)
%
%   INPUTS
%        H,k: halfspace representation variables polyhedron
%
%   OUTPUT
%    Hrm,krm: the resulting polyhedron

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

H_ = repmat({H},N,1);
Hrm = blkdiag(H_{:});
krm = repmat(k,N,1);

end
