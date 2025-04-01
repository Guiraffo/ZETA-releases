function ret = abs(a)
%abs(a) returns the absolute value of interval a

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
    
if(numel(a)>1)
    error('This function is only implemented for scalars.');
end    
    
aL = inf(a);
aU = sup(a);

if(aL>=0)
    ret = a;
elseif(aU<=0)
    ret = -a;
else
    ret = Interval(0,max([abs(a.LB),abs(a.UB)]));
end
    
end