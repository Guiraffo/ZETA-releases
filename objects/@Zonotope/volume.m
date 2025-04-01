function Vol = volume(Z,method)
%VOLUME computes the volume of a zonotope
%
%   SYNTAX: Vol = VOLUME(Z)
%           Vol = VOLUME(Z,method)
%
%   INPUTS
%           Z: zonotope object
%      method: 'Hrep': converts Z to H-rep and compute the volume using MPT
%              'Grep': computes the volume using Eq.(6) in Alamo et al
%                      (2005)
%           'Grepvpa': computes the volume using Eq.(6) in Alamo et al
%                      (2005), with variable precision arithmetic (slow)
%           'partope': computes the approximated volume of Z, by
%                      calculating the volume in G-rep of the parallelotope
%                      obtained after using generator reduction
%               'box': computes the volume of the interval hull of Z
%              (default is 'Hrep')
%
%   OUTPUT
%         Vol: the (approximated) volume of Z computed using the chosen
%              method

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

Z = Zonotope(Z); % Redundant, but for safety

switch method
    case 'Hrep'
        
        Vol = volume(Polyhedron(Z));
        
    case 'Grep'
       
        % Requires the 'onecomb' matlab function available online

        dim_z = size(Z.G,1);
        nof_g = size(Z.G,2);
        %Combs = nchoosek(1:nof_g,dim_z);

        %nof_combs = size(Combs,1);
        nof_combs = nchoosek(nof_g,dim_z);

        if(nof_combs > 1.0e6)
            disp(['Warning: ', num2str(nof_combs), ' possible combinations in Zonotope/volume'])
        end

        Vol = 0;
        for i=1:nof_combs
            Vol = Vol + abs(det(Z.G(:,thirdparty.onecomb(nof_g,dim_z,i))));
        end
        Vol = (2^dim_z)*Vol;      
      
    case 'Grepvpa'
        
        dim_z = size(Z.G,1);
        nof_g = size(Z.G,2);

        nof_combs = nchoosek(nof_g,dim_z);

        if(nof_combs > 1.0e6)
            disp(['Warning: ', num2str(nof_combs), ' possible combinations in volume'])
        end

        Vol = 0;
        for i=1:nof_combs
            Vol = Vol + abs(det(vpa(Z.G(:,thirdparty.onecomb(nof_g,dim_z,i)))));
        end
        Vol = (2^dim_z)*Vol;        
        
        
    case 'partope'
        
        Vol = volume(reduction(Z,size(Z.c,1)),'Grep');
        
    case 'box'
        
        Vol = prod(diam(intervalhull(Z)));        
        
    otherwise
        error('Invalid method in Zonotope/volume.')
end

end

