function Z = removenull(Zin)
%REMOVENULL removes null (zero) segments from a zonotope
%
%   SYNTAX: Z = REMOVENULL(Zin)
%
%   INPUTS
%         Zin: zonotope object
%
%   OUTPUT
%           Z: zonotope without null segments

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

r = size(Zin.G,2);
H = Zin.G;

i = 1;
while i<=r

    if(norm(H(:,i),1) == 0)
        
        H(:,i) = [];
        r = r - 1;
    
    else
        i = i + 1;
    end
    
end

Z = Zonotope(Zin.c, H);

end

