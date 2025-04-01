function Z = abs(X)
%ABS returns the absolute value a Polyrelax object
%
%   SYNTAX: Z = ABS(X)
%
%   INPUTS
%           Z: the Polyrelax object
%
%   OUTPUT
%           Z: the resulting Polyrelax object

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

xL = inf(X.x);
xU = sup(X.x);

% Trivial cases
if(xU<=0)
    Z = -X;
    return;
elseif(xL>=0)
    Z = X;
    return;
end

% If its not a trivial case (i.e., input contains zero), do

Z = Polyrelax(abs(X.x)); 

% Initialize polyhedral enclosure
Hadd = zeros(3,Z.i);
kadd = zeros(3,1);

% Convex lowerbound
% z > x;
% z > -x;
Hadd(1,X.i) = 1;
Hadd(2,X.i) = -1;
Hadd(1,Z.i) = -1;
Hadd(2,Z.i) = -1;

% Concave upperbound (secant)
% Obtained from z < (abs(xU) - abs(xL)/(xU-xL) + abs(xL), with xU positive and
% xL is negative
% z < (xU + xL)/(xU - xL) - xL
concv_slope = (xU + xL)/(xU - xL);
Hadd(3,X.i) = -concv_slope;
Hadd(3,Z.i) = 1; 
kadd(3,1) = -concv_slope*xL - xL;

% Adds the inequality constraints to Polyrelax Hrep
Polyrelax.addHrepHrow(Hadd,kadd);        
    
end
