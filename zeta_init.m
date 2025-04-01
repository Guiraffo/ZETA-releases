function zeta_init
% ZETA_INIT initializes ZETA

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

clear all;

disp(' ');
disp('============================== <strong>ZETA</strong> ================================');
disp('========== <strong>Z</strong>onotope-based <strong>E</strong>s<strong>T</strong>imation and f<strong>A</strong>ult diagnosis ===========');
disp('==================== of discrete-time systems ======================');
disp(' ');
disp('Checking for 3rd-party libraries...');
disp(' ');

% Adds folder and subfolders to MATLAB path
addpath(erase(mfilename('fullpath'), mfilename));
addpath([erase(mfilename('fullpath'), mfilename),'objects']);
addpath([erase(mfilename('fullpath'), mfilename),'packages']);
addpath([erase(mfilename('fullpath'), mfilename),'estimation']);
addpath([erase(mfilename('fullpath'), mfilename),'faultdiag']);

% Checks for Yalmip
disp('[?] Trying to find Yalmip.');
if(exist('yalmipdemo','file')~=2)
    error('[X] Could not find Yalmip. Yalmip is required to run ZETA.');
else
    disp('[+] Found Yalmip.');
end
% Checks for Gurobi 
disp('[?] Trying to find Gurobi.');
if(exist('gurobi','file')~=3)
    disp('[X] Could not find Gurobi. Most of the optimization problems will still be solved using MATLAB Optimization Toolbox.');
    toolsettings.solverlp('linprog');
    toolsettings.solvermilp('intlinprog');
else
    disp('[+] Found Gurobi.');
    toolsettings.solverlp('gurobi');
    toolsettings.solvermilp('gurobi');    
end
% Checks for MPT
disp('[?] Trying to find MPT.');
if(exist('mpt_init','file')~=2)
    disp('[X] Could not find MPT. All set plotting will be handled by Yalmip, but some operations may be unavailable.');
    toolsettings.plot('yalmip');
else
    disp('[+] Found MPT.');
end

disp(' ');
disp('Checking finished!');

toolsettings.show;




end