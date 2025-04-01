function [lambda,min_eigenv,Hessian_vert] = getlambda_alphaBB(f,Zhull,varargin)
% GETLAMBDA_ALPHABB implements the method in eq. (12) by Adjiman and
%                   Floudas (1996) to obtain the lambda parameter of the
%                   convex underestimator. For internal use in DC
%                   programming propagations only.

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

% Check inputs
switch nargin
    case 2 % Called in propagation
        Hessian = f.Hessian(Zhull);
    case 3 % Called in state estimation, update
        nof_x = varargin{1};
        Hessian = f.Hessian(Zhull(1:nof_x),Zhull(nof_x+1:end));
    case 5 % Called in state estimation, prediction
        nof_x = varargin{1};
        u = varargin{2};
        Ts = varargin{3};
        Hessian = f.Hessian(Zhull(1:nof_x),Zhull(nof_x+1:end),u,Ts);       
    otherwise
        error('Invalid call of getlambda_alphaBB.');
end

nof_z = size(Zhull,1);        
nof_f = length(Hessian);

% Get the first half set of vertices of the unitary box
nof_vert = 2^(nof_z-1);
UnitB_vert = util.box_vertices(Interval(-ones(nof_z,1),ones(nof_z,1)),nof_vert); 

% Rows are the rows of f, columns are the vertices
Hessian_vert = cell(nof_f,nof_vert);
min_eigenv = zeros(nof_f,nof_vert);

% Compute the 2^(nz-1) vertices of the hessian matrices and their lowest eigenvalues
for m=1:nof_f
    Hessian{m} = Interval(Hessian{m}); % Interval call in case its not an interval matrix
    for k=1:nof_vert
        Hessian_vert{m,k} = zeros(nof_z,nof_z);
        for i=1:nof_z
            for j=1:nof_z
                if(i==j)
                    Hessian_vert{m,k}(i,j) = inf(Hessian{m}(i,j));
                else
                    signprod = UnitB_vert(i,k)*UnitB_vert(j,k);
                    if(signprod>=0)
                        Hessian_vert{m,k}(i,j) = inf(Hessian{m}(i,j));
                    elseif(signprod<0)
                        Hessian_vert{m,k}(i,j) = sup(Hessian{m}(i,j));
                    end
                end
            end
        end
        eigenv = eig(Hessian_vert{m,k});
        min_eigenv(m,k) = min(eigenv);
    end
end

% Minimum eigenvalues of the interval hessian matrices
min_eigenv = min(min_eigenv,[],2);

lambda = max([zeros(nof_f,1),-min_eigenv],[],2);
    
end

