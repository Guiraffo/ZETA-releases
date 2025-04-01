function S = intersection(S1,S2)
%INTERSECTION computes the intersection of two parallel strips
%
%   SYNTAX: S = INTERSECTION(S1,S2)
%
%   INPUTS
%           S1: first strip object
%           S2: second strip object
%
%   OUTPUT
%           S: strip corresponding to the intersection

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

if(S1.p == S2.p)
    
    switch toolsettings.IA_class
        case 'intval'
        
            Interval1 = midrad(S1.d, S1.sigma);
            Interval2 = midrad(S2.d, S2.sigma);

            Intersection = intersect(Interval1, Interval2);

            d = mid(Intersection);
            sigma = rad(Intersection);

            if(isnan(d)||isnan(sigma)) % If the strips do not intersect each other

                %error('Error in strip_intersection: Strips do not intersect each other!');
                error('Empty intersection detected!');

            else
                    S = Strip(S1.p, d, sigma);
            end
            
        case 'Interval'
            Interval1 = Interval(S1.d - S1.sigma, S1.d + S1.sigma);
            Interval2 = Interval(S2.d - S2.sigma, S2.d + S2.sigma);
            [Intersection,emptyint] = intersect(Interval1,Interval2);
            if(emptyint)
                error('Empty intersection detected!');
            else
                d = mid(Intersection);
                sigma = rad(Intersection);                
                S = Strip(S1.p, d, sigma);
            end
        otherwise
            error('Invalid IA class.');
    end
            
    
else
    
    error('Strips are not parallel!');
    
end



end

