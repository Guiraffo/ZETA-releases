function [ret,isempty] = intersect(a,b)
%INTERSECT returns the intersection of two intervals

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

a = Interval(a);
b = Interval(b);

dima1 = size(a,1);
dima2 = size(a,2);
dimb1 = size(b,1);
dimb2 = size(b,2);

if((dima1~=dimb1)||(dima2~=dimb2))
    error('Input dimensions mismatch.');
end

aL = inf(a);
aU = sup(a);
bL = inf(b);
bU = sup(b);

ret = Interval.zeros(dima1,dima2);

if(nargout==2)
    isempty = 0;
end    
for i=1:dima1
    for j=1:dima2
        if ((aL(i,j) <= bU(i,j)) && (bL(i,j) <= aU(i,j)))
            ret(i,j) = Interval(max([aL(i,j),bL(i,j)]), min([aU(i,j),bU(i,j)]));
        else
            if(nargout==2)
                ret = Interval.empty;
                isempty = 1;        
            else
                error('Intersection is an empty set.');
            end
        end
    end
end
    
end