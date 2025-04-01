function plot(Z,color,alpha,borderstyle,canvas)
%PLOT plots a line zonotope bounded by either a given canvas or a multiple
%     of the current figure axis settings
%
%   SYNTAX: LZONOTOPE(Z,color,alpha)
%           LZONOTOPE(Z,color,alpha,borderstyle)
%           LZONOTOPE(Z,color,alpha,borderstyle,canvas)
%
%   INPUTS
%             Z: mixed zonotope as structure object
%         color: color with which the line zonotope will be filled
%                (RGB vector or a string name of the color, e.g. 'blue')
%         alpha: opacity of the line zonotope (between 0 and 1)
%   borderstyle: line style of the set borders (default is '-')
%        canvas: canvas limits for plotting unbounded sets (same format as
%                axis function)

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
    canvas = axis;
elseif(nargin==2)
    alpha = 1;
    borderstyle = '-';      
    canvas = axis;    
elseif(nargin==3)
    borderstyle = '-';  
    canvas = axis;
elseif(nargin==4)
    canvas = axis;    
end

% Scale canvas
canvasdim = length(canvas)/2;
canvas = reshape(canvas,2,canvasdim).';
canvasmid = 0.5*(canvas(:,2) + canvas(:,1));
canvasrad = 0.5*(canvas(:,2) - canvas(:,1));
canvas(:,1) = canvasmid - 3*canvasrad;
canvas(:,2) = canvasmid + 3*canvasrad;


lzonotope_dimension = Z.dim;
if(lzonotope_dimension > 3)
    error('This function plots line zonotopes only up to dimension three.');
end

switch toolsettings.plot
    case {'yalmip','mpt-H','mpt-V'}
        
        switch lzonotope_dimension
            case 1
                nof_points = 10;
            case 2
                nof_points = 500;
            case 3
                nof_points = 2000;
            case 0
                error('Error in the lzonotope representation LZonotope/plot.');
            otherwise
                error('This function plots line zonotopes only up to dimension three.');
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
                error('Invalid color string in LZonotope/plot.')
        end
        
        nof_generators = size(Z.G,2);
        nof_ugenerators = size(Z.M,2);

        x = sdpvar(lzonotope_dimension,1);
        xi = sdpvar(nof_generators,1);
        delta = sdpvar(nof_ugenerators,1);

        OPTIONS = sdpsettings;
        OPTIONS.plot.lighting = 0;
        OPTIONS.plot.shade = alpha;
        OPTIONS.plot.waitbar = 0;
        OPTIONS.plot.wirestyle = borderstyle;        
        plot([x == Z.c + Z.G*xi + Z.M*delta; norm(xi,Inf) <= 1; Z.S*delta + Z.A*xi == Z.b; canvas(:,1) <= x; x <= canvas(:,2)],x,color,nof_points,OPTIONS);
        grid on;          

    otherwise
        error('Invalid method in LZonotope/plot.')
end

end
