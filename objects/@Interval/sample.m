function Sample = sample(B,N,method)
%SAMPLE generate sample points enclosed by an interval
%
%   SYNTAX: Sample = SAMPLE(B,N)
%           Sample = SAMPLE(B,N,method)
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

if(nargin<2)
    error('Not enough input arguments.');
elseif(nargin==2)
    method = 'uniform';
end

dim  = size(B,1);
dim2 = size(B,2);
if(dim2~=1)
    error('Interval/sample only works for column interval vectors.')
end


switch method
    case 'uniform' % Uniform sampling
        
        Sample = repmat(mid(B),1,N) + diag(rad(B))*(2*rand(dim,N) - ones(dim,N));
    
    case 'grid' % Uniform grid
            
        % Gridding values
        Grid = cell(1,dim);
        for j=1:dim
            Grid{j} = linspace(inf(B(j)),sup(B(j)),N+1);
        end

        % Grid samples
        Sample_coord = cell(1,dim);
        [Sample_coord{:}] = ndgrid(Grid{:});

        Sample = zeros(dim, (N+1)^dim);
        for j=1:dim
            Sample(j,:) = reshape(Sample_coord{j},1,(N+1)^dim);
        end  
        
    case 'Gaussian' % Gaussian sampling     
        
        Sample = zeros(dim,N);
        
        % Sampling by rejection
        j = 1;
        while(j<=N)

            % Sampling in the interval hull
            gausssample = mid(B) + diag(rad(B)/3)*randn(dim,1);

            if(isinside(B,gausssample))
                Sample(:,j) = gausssample;
                j = j + 1;
            end    

        end     
        
    otherwise
        error('Invalid method in Interval/sample.');
end


end

