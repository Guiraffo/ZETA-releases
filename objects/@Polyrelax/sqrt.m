function Z = sqrt(X)
%SQRT returns the square root a polyrelax object
%
%   SYNTAX: Z = SQRT(X)
%
%   INPUTS
%           Z: the polyrelax object
%
%   OUTPUT
%           Z: the resulting polyrelax object

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

if(~(xL>=0))
    error('Only positive inputs are allowed.')
end

Z = Polyrelax(sqrt(X.x));

% Tolerance value for 'near zero' lower bound
TOL = 1e-2;

% Compute the envelopes

xmid = 0.5*(xL+xU);
Hadd = zeros(4,Z.i);
kadd = zeros(4,1);

% Convex lowerbound
% z > ((sqrt(xU) - sqrt(xL))/(xU - xL))*(x - xL) + sqrt(xL)
convx_slope = ((sqrt(xU) - sqrt(xL))/(xU - xL));
Hadd(1,X.i) = convx_slope;
Hadd(1,Z.i) = -1; 
kadd(1,1) = convx_slope*xL - sqrt(xL);

% Concave upperbound
% IF xL > TOL, {z < diff(sqrt(x),x)|(x=xL)*(x - xL) + sqrt(xL)}, ELSE, {x > xL} (sqrt(x) does not have a derivative at x=0)
% z < diff(sqrt(x),x)|(x=xmid)*(x - xmid) + sqrt(xmid)
% z < diff(sqrt(x),x)|(x=xU)*(x - xU) + sqrt(xU)
if(xL>TOL)
    Hadd(2,X.i) = -diffsqrt(xL);
    Hadd(2,Z.i) = 1;
    kadd(2,1) = -diffsqrt(xL)*xL + sqrt(xL);
else
    Hadd(2,X.i) = -1;
    kadd(2,1) = -xL;
end
Hadd(3,X.i) = -diffsqrt(xmid);
Hadd(4,X.i) = -diffsqrt(xU);
Hadd(3,Z.i) = 1;
Hadd(4,Z.i) = 1;
kadd(3,1) = -diffsqrt(xmid)*xmid + sqrt(xmid);
kadd(4,1) = -diffsqrt(xU)*xU + sqrt(xU);




    
% Adds the inequality constraints to polyrelax Hrep
Polyrelax.addHrepHrow(Hadd,kadd);

end

function out = diffsqrt(x)
out = 0.5/sqrt(x);
end