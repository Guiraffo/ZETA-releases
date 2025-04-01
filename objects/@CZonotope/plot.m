function plot(Z,color,alpha,borderstyle)
%PLOT plots a constrained zonotope
%
%   SYNTAX: PLOT(Z)
%           PLOT(Z,color)
%           PLOT(Z,color,alpha)
%           PLOT(Z,color,alpha,borderstyle)
%
%   INPUTS
%             Z: constrained zonotope object
%         color: color with which the set will be filled (RGB vector or a
%                string name of the color, e.g. 'blue') (default is 'red')
%         alpha: opacity of the set (between 0 and 1) (default is 1)
%   borderstyle: line style of the set borders (default is '-')

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

if(nargin==0)
    error('Not enough input arguments.');
elseif(nargin==1)
    color = 'red';
    alpha = 1;
    borderstyle = '-';      
elseif(nargin==2)
    alpha = 1;
    borderstyle = '-';      
elseif(nargin==3)
    borderstyle = '-';  
end


czonotope_dimension = Z.dim;
if(czonotope_dimension > 3)
    error('This function plots contrained zonotopes only up to dimension three.');
end

switch toolsettings.plot
    case {'mpt-H','mpt-V'}
        plot(Polyhedron(Z),'Color',color,'Alpha',alpha,'LineStyle',borderstyle);          
    case 'yalmip'
        
        switch czonotope_dimension
            case 1
                nof_points = 10;
            case 2
                nof_points = 500;
            case 3
                nof_points = 2000;
            case 0
                error('Error in CZonotope/plot: the input constrained zonotope has dimension zero.');
            otherwise
                error('This function plots constrained zonotopes only up to dimension three.');
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
                error('Invalid color string in CZonotope/plot.')
        end
        
        nof_generators = size(Z.G,2);

        x = sdpvar(czonotope_dimension,1);
        xi = sdpvar(nof_generators,1);

        OPTIONS = sdpsettings;
        OPTIONS.plot.lighting = 0;
        OPTIONS.plot.shade = alpha;
        OPTIONS.plot.waitbar = 0;
        OPTIONS.plot.wirestyle = borderstyle;
        plot([x == Z.c + Z.G*xi; norm(xi,Inf) <= 1; Z.A*xi == Z.b],x,color,nof_points,OPTIONS);
        grid on;   
              
    otherwise
        error('Invalid method in CZonotope/plot.')
end

end
