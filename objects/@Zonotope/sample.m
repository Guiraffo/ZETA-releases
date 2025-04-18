function Sample = sample(Z,N,method)
%SAMPLE generate sample points enclosed by a zonotope
%
%   SYNTAX: Sample = SAMPLE(Z,N)
%           Sample = SAMPLE(Z,N,method)
%
%   INPUTS
%           Z: zonotope object
%           N: number of samples
%      method: 'uniform': generates N uniformly distributed samples in Z
%                 'grid': generates a NxN grid with samples in Z (slow)
%             'Gaussian': generates N samples in Z according to a truncated
%                         Gaussian distribution with mean equal to the
%                         center of Z, and standard deviation equal to 1/3
%                         of the component-wise radius of Z
%                         (default method is 'uniform')
%
%   OUTPUT
%      Sample: generated samples

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

dim = size(Z.c,1);

if(nargin<2)
    error('Not enough input arguments.');
elseif(nargin==2)
    method = 'uniform';
end


switch method
    case 'uniform' % Uniform sampling
        
        Sample = zeros(dim,N);
        Hull = intervalhull(Z);

        % Sampling by rejection
        j = 1;
        while(j<=N)


            % Sampling in the interval hull
            boxsample = mid(Hull) + diag(rad(Hull))*(2*rand(dim,1) - ones(dim,1));

            if(isinside(Z,boxsample))
                Sample(:,j) = boxsample;
                j = j + 1;
            end    

        end
    
    case 'grid' % Uniform grid
        
        Hull = intervalhull(Z);
        
        % Gridding values
        Grid = cell(1,dim);
        for j=1:dim
            %Grid(j,:) = linspace(inf(Hull(j)),sup(Hull(j)),N+1);
            Grid{j} = linspace(inf(Hull(j)),sup(Hull(j)),N+1);
        end

        % Grid samples
        Sample_coord = cell(1,dim);
        [Sample_coord{:}] = ndgrid(Grid{:});

        Sample = zeros(dim, (N+1)^dim);
        for j=1:dim
            Sample(j,:) = reshape(Sample_coord{j},1,(N+1)^dim);
        end

        % Remove samples that do not belong to the zonotope
        i = 1;
        nof_samples = size(Sample,2);
        while(i<=nof_samples)   
            if(~isinside(Z,Sample(:,i)))
                Sample(:,i) = [];
                nof_samples = nof_samples - 1;
            else
                i = i + 1;
            end  
        end   
        
    case 'Gaussian' % Gaussian sampling     
        
        Sample = zeros(dim,N);
        Hull = intervalhull(Z);

        % Sampling by rejection
        j = 1;
        while(j<=N)

            % Sampling in the interval hull
            boxsample = mid(Hull) + diag(rad(Hull)/3)*randn(dim,1);

            if(isinside(Z,boxsample))
                Sample(:,j) = boxsample;
                j = j + 1;
            end    

        end     
        
    otherwise
        error('Invalid method in Zonotope/sample.');
end

end

