function disp(Z)
%DISP prints a line zonotope in the command window
%
%   SYNTAX: DISP(Z)
%
%   INPUTS
%           Z: line zonotope

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

disp('Line zonotope with properties:');
disp('c =');   if(isempty(Z.c)); disp(['empty(',num2str(size(Z.c,1)),',',num2str(size(Z.c,2)),')']); else; disp(Z.c); end;
disp('G =');   if(isempty(Z.G)); disp(['empty(',num2str(size(Z.G,1)),',',num2str(size(Z.G,2)),')']); else; disp(Z.G); end;
disp('M =');   if(isempty(Z.M)); disp(['empty(',num2str(size(Z.M,1)),',',num2str(size(Z.M,2)),')']); else; disp(Z.M); end;
disp('S =');   if(isempty(Z.S)); disp(['empty(',num2str(size(Z.S,1)),',',num2str(size(Z.S,2)),')']); else; disp(Z.S); end;
disp('A =');   if(isempty(Z.A)); disp(['empty(',num2str(size(Z.A,1)),',',num2str(size(Z.A,2)),')']); else; disp(Z.A); end;
disp('b =');   if(isempty(Z.b)); disp(['empty(',num2str(size(Z.b,1)),',',num2str(size(Z.b,2)),')']); else; disp(Z.b); end;

end



