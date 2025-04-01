function [Z,P] = inclusion(X,J,Xbar)
%INCLUSION computes the CZ-inclusion: a constrained zonotope that contains
%          a set defined as the product of an interval matrix and a
%          constrained zonotope, S = JX (Theorem 1 in Rego et al. (2020))
%
%   SYNTAX: Z = INCLUSION(X,J)
%           Z = INCLUSION(X,J,Xbar)
%       [Z,P] = INCLUSION(X,J)
%       [Z,P] = INCLUSION(X,J,Xbar)
%
%   INPUTS
%           X: constrained zonotope object
%           J: interval matrix
%        Xbar: zonotope Xbar to be used in the CZ-inclusion
%              (optional. If not given, Xbar will be computed by 
%               constraint elimination)
%
%   OUTPUT
%           Z: constrained zonotope containing S = JX
%           P: generator matrix of the box part

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

dimension_s = size(X.c,1); % space dimension
n_g = size(X.G,2); % number of generators
n_c = size(X.A,1); % number of constraints

if(nargin<2)
    error('Not enough input arguments.');
elseif(nargin==2)  
    % If the zonotope Xbar is must be computed  
    
    % If the n > n_g-n_c, throws an error (the obtained zonotope is not full dimensional)
    if((dimension_s > n_g-n_c))
        error('The obtained zonotope by constraint elimination is not full dimensional!');
    end

    % Builds a cell array for linking with libCZon
    X_{1} = X.c;
    X_{2} = X.G;
    X_{3} = X.A;
    X_{4} = X.b;

    % Obtain Xbar by eliminating all the constraints in X
    for j=1:n_c
        X_ = libCZon.ScaleDualize(X_);
        
        % Stop condition due to the possibility of removing null constraints inside libCZon.ScaleDualize
        if(size(X_{4},1)==0)
            break;
        end        
        
    end

    pbar = X_{1};
    Mbar = X_{2};

    m = (J - mid(J))*pbar;
    midm = mid(m);
    diamm = diam(m);

    p = X.c;
    M = X.G;
    A = X.A;
    b = X.b;

    Gbar = (1/2)*diag(sum(diam(J)*abs(Mbar),2));

    P = (1/2)*diag(diamm) + Gbar;

    Z = CZonotope(midm + mid(J)*p, [mid(J)*M, P], [A, zeros(size(A,1),size(J,1))], b);
    
else
    
    % Xbar is given
    
    pbar = Xbar.c;
    Mbar = Xbar.G;

    m = (J - mid(J))*pbar;
    midm = mid(m);
    diamm = diam(m);

    p = X.c;
    M = X.G;
    A = X.A;
    b = X.b;

    Gbar = (1/2)*diag(sum(diam(J)*abs(Mbar),2));

    P = (1/2)*diag(diamm) + Gbar;

    Z = CZonotope(midm + mid(J)*p, [mid(J)*M, P], [A, zeros(size(A,1),size(J,1))], b);

end
    
   
end

