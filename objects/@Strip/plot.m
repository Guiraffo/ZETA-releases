function plot(S,color,alpha,borderstyle,canvas)
%PLOT plots a strip
%
%   SYNTAX: PLOT(S,color,alpha)
%   SYNTAX: PLOT(S,color,alpha,borderstyle)
%   SYNTAX: PLOT(S,color,alpha,borderstyle,canvas)
%
%   INPUTS
%           S: strip in struct format
%       color: color with which the strip will be filled (RGB vector or a
%              string name of the color, e.g. 'blue')
%       alpha: opacity of the strip (between 0 and 1)
% borderstyle: line style of the set borders (default is '-')
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

strip_dimension = size(S.p,1);
if(strip_dimension > 3)
    error('This function plots strip only up to dimension three.');
end

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

% Save current axis
axis0 = axis;  

switch toolsettings.plot
    case 'yalmip'

        switch strip_dimension
            case 1
                nof_points = 10;
            case 2
                nof_points = 500;
            case 3
                nof_points = 2000;
            case 0
                error('Error in the strip representation Strip/plot.');
            otherwise
                error('This function plots strips only up to dimension three.');
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
                error('Invalid color string in Strip/plot.')
        end
        
        x = sdpvar(strip_dimension,1);

        OPTIONS = sdpsettings;
        OPTIONS.plot.lighting = 0;
        OPTIONS.plot.shade = alpha;
        OPTIONS.plot.waitbar = 0;
        OPTIONS.plot.wirestyle = borderstyle;        
        plot([S.p'*x - S.d >= -S.sigma; S.p'*x - S.d <= S.sigma;  canvas(:,1) <= x; x <= canvas(:,2)],x,color,nof_points,OPTIONS);
        grid on;              
        
    case {'mpt-H','mpt-V'}
        CanvasBox = Polyhedron('LB',canvas(:,1),'UB',canvas(:,2));
        plot(intersect(Polyhedron(S),CanvasBox),'Color',color,'Alpha',alpha,'LineStyle',borderstyle);    
    otherwise
        error('Invalid method in Strip/plot.');    
end

axis(axis0);

end
