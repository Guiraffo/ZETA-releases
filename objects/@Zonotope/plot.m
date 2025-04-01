function plot(Z,color,alpha,borderstyle)
%PLOT plots a zonotope
%
%   SYNTAX: PLOT(Z)
%           PLOT(Z,color)
%           PLOT(Z,color,alpha)
%           PLOT(Z,color,alpha,borderstyle)
%
%   INPUTS
%             Z: zonotope object
%         color: color with which the set will be filled (RGB vector or
%                a string name of the color, e.g. 'blue') (default is 'red')
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


zonotope_dimension = Z.dim;
if(zonotope_dimension > 3)
    error('This function plots zonotopes only up to dimension three.');
end

switch toolsettings.plot
    case 'mpt-H'     
        plot(Polyhedron(Z,method),'Color',color,'Alpha',alpha,'LineStyle',borderstyle);       
    case 'mpt-V'    
        plot(Polyhedron(Z,method),'Color',color,'Alpha',alpha,'LineStyle',borderstyle);   
    case 'yalmip'
        
        switch zonotope_dimension
            case 1
                nof_points = 10;
            case 2
                nof_points = 500;
            case 3
                nof_points = 2000;
            case 0
                error('Error in zonotope/plot: the input zonotope has dimension zero.');
            otherwise
                error('This function plots zonotopes only up to dimension three.');
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
                error('Invalid color string in zonotope/plot.')
        end
        
        nof_generators = size(Z.G,2);

        x = sdpvar(zonotope_dimension,1);
        xi = sdpvar(nof_generators,1);

        OPTIONS = sdpsettings;
        OPTIONS.plot.lighting = 0;
        OPTIONS.plot.shade = alpha;
        OPTIONS.plot.waitbar = 0;
        OPTIONS.plot.wirestyle = borderstyle;
        plot([x == Z.c + Z.G*xi; norm(xi,Inf) <= 1],x,color,nof_points,OPTIONS);
        grid on; 
        
        
    otherwise
        error('Invalid method in zonotope/plot.')
end


% Revision 17-05-20: merged with plot_szonotope_v2. Added also the 'Yalmip' option.

end
