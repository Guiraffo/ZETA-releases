function f = fhandle(fname)
%FHANDLE generates function handles for a given function name, necessary
%        for propagation of sets through nonlinear functions and nonlinear
%        state estimation. Internal use only.

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

% Common function handles
f.eval = str2func(fname);
f.Jacob = str2func(strcat([fname,'_Jacob']));
f.HessianTriang = str2func(strcat([fname,'_HessianTriang']));
f.Hessian = str2func(strcat([fname,'_Hessian']));

% DC programming function handles for exact decomposition
f.dcA = str2func(strcat([fname,'_dcA']));
f.dcA_Jacob = str2func(strcat([fname,'_dcA_Jacob']));    
f.dcB = str2func(strcat([fname,'_dcB']));    
f.dcB_Jacob = str2func(strcat([fname,'_dcB_Jacob']));        


end    