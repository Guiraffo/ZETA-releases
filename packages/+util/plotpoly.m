function plotpoly(H,k,color,alpha,borderstyle)
%PLOTPOLY plots a polyhedron in halfspace representation
%
%   SYNTAX: PLOTPOLY(S,color,alpha)
%           PLOTPOLY(S,color,alpha,borderstyle)
%
%   INPUTS
%           S: strip in struct format
%       color: color with which the strip will be filled (RGB vector or a
%              string name of the color, e.g. 'blue')
%       alpha: opacity of the strip (between 0 and 1)
% borderstyle: line style of the set borders (default is '-')

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

poly_dimension = size(H,2);
if(poly_dimension > 3)
    error('This function plots polyhedrons only up to dimension three.');
end

if(nargin<2)
    error('Not enough input arguments.');
elseif(nargin==2)
    color = 'red';
    alpha = 1;
    borderstyle = '-';
elseif(nargin==3)
    alpha = 1;
    borderstyle = '-';          
elseif(nargin==4)
    borderstyle = '-';   
end


switch toolsettings.plot
    case 'yalmip'

        switch poly_dimension
            case 1
                nof_points = 10;
            case 2
                nof_points = 500;
            case 3
                nof_points = 2000;
            case 0
                error('Error in the polyhedron representation in plotpoly.');
            otherwise
                error('This function plots polyhedrons only up to dimension three.');
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
                error('Invalid color string in plotpoly.')
        end
        
        x = sdpvar(poly_dimension,1);

        OPTIONS = sdpsettings;
        OPTIONS.plot.lighting = 0;
        OPTIONS.plot.shade = alpha;
        OPTIONS.plot.waitbar = 0;
        OPTIONS.plot.wirestyle = borderstyle;        
        plot([H*x <= k],x,color,nof_points,OPTIONS);
        grid on;              
        
    case {'mpt-H','mpt-V'}
        plot(Polyhedron(H,k),'Color',color,'Alpha',alpha,'LineStyle',borderstyle);    
    otherwise
        error('Invalid method in plotpoly.');    
end


end
