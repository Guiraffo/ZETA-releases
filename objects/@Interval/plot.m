function plot(B,color,alpha)
%PLOT plots an interval
%
%   SYNTAX: PLOT(B,color,alpha)
%
%   INPUTS
%           B: box as an interval column vector
%       color: color with which the box will be filled (RGB vector or a
%              string name of the color, e.g. 'blue')
%       alpha: opacity of the box (between 0 and 1)

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

box_dimension = size(B,1);
if(box_dimension > 3)
    error('This function plots interval vectors only up to dimension three.');
end
if(size(B,2)>1)
    error('This function plots only column interval vectors.');
end

switch toolsettings.plot
    case {'mpt-H','mpt-V'}
        plot(Polyhedron(B),'Color',color,'Alpha',alpha);      
    case 'yalmip'
        switch box_dimension
            case 1
                nof_points = 10;
            case 2
                nof_points = 500;
            case 3
                nof_points = 2000;
            case 0
                error('Error in interval/plot: the input interval has dimension zero.');
            otherwise
                error('This function plots interval vectors only up to dimension three.');
        end

        % Color string wrapper to yalmip plot (this is different from MPT plot)
        switch color
            case 'red'
                color = 'r';
            case 'yellow'
                color = 'y';
            case 'magenta'
                color = 'm';
            case 'cyan'
                color = 'c';
            case 'green'
                color = 'g';                
            case 'blue'
                color = 'b';
            case 'black'
                color = 'k';
            otherwise
                error('Invalid color string in interval/plot.')
        end
        
        
        x = sdpvar(box_dimension,1);
        LB = inf(B);
        UB = sup(B);
        
        OPTIONS = sdpsettings;
        OPTIONS.plot.lighting = 0;
        OPTIONS.plot.shade = alpha;
        OPTIONS.plot.waitbar = 0;
        if(nargin==5)
            OPTIONS.plot.wirestyle = borderstyle;
        end
        plot([LB <= x; x <= UB],x,color,nof_points,OPTIONS);
        grid on;         
        
    otherwise
        error('Invalid plot method.');
end

end
