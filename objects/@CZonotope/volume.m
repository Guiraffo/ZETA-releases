function Vol = volume(Z,method)
%VOLUME computes the volume of a constrained zonotope
%
%   SYNTAX: Vol = VOLUME(Z)
%           Vol = VOLUME(Z,method)
%
%   INPUTS
%           Z: constrained zonotope object
%      method: 'Hrep': converts Z to H-rep and compute the volume (requires
%                      MPT)
%           'partope': computes the volume of a paralelotope bound of Z
%        'partopevpa': computes the volume of a paralelotope bound of Z,
%                      using variable precision arithmetic (slow)
%   'partope-nthroot': computes the nthroot of the volume of a paralelotope
%                      bound of Z
%               'box': computes the volume of the interval hull of Z
%       'box-nthroot': computes the nthroot of the volume of a interval
%                      hull of Z
%              (default is 'Hrep')
%
%   OUTPUT
%         Vol: the resulting volume

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

if(nargin<1)
    error('Not enough input arguments.');
elseif(nargin==1)
    method = 'Hrep';
end

% Check for empty input
if(isempty(Z))
    Vol = 0;
    return;
end

% Check for degenerate set
if((size(Z.G,2)-size(Z.A,1))<size(Z.c,1))
    Vol = 0;
    return;
end


switch method
    case 'Hrep'
        
        Vol = volume(Polyhedron(Z));
                  
    case 'partope'
        
        Zpartope = partopebound(Z,'rescale');
        Vol = volume(Zonotope(Zpartope.c,Zpartope.G),'Grep');
        
    case 'partopevpa'
        
        Zpartope = partopebound(Z,'rescale');
        Vol = volume(Zonotope(Zpartope.c,Zpartope.G),'Grepvpa');        
        
    case 'partope-nthroot'
        
        Zpartope = partopebound(Z,'rescale');
        Vol = nthroot(volume(Zonotope(Zpartope.c,Zpartope.G),'Grep'),Z.dim);
        
    case 'box'
        
        Vol = prod(diam(intervalhull(Z)));
        
    case 'box-nthroot'
        
        Vol = nthroot(prod(diam(intervalhull(Z))),Z.dim);
        
    otherwise
        error('Invalid method in CZonotope/volume.')
end
        
end

