function Z = sec(X)
%SEC returns the secant of a polyrelax object
%
%   SYNTAX: Z = SEC(X)
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

% Input check
xmin = -pi/2;
xmax = pi/2;
xL = inf(X.x);
xU = sup(X.x);

% Check if input is not contained in [-pi/2,pi/2]
if((xL<=xmin)||(xL>=xmax)||(xU<=xmin)||(xU>=xmax))
    error(' The secant function is currently implemented only for the interval [-pi/2,pi/2].');
end

Z = Polyrelax(sec(X.x)); 
   

% The secant function is convex in [-pi/2,pi/2]
        
xmid = 0.5*(xL+xU);
Hadd = zeros(4,Z.i);
kadd = zeros(4,1);

% Concave upperbound
% z < ((sec(xU) - sec(xL))/(xU - xL))*(x - xL) + sec(xL)
concv_slope = ((sec(xU) - sec(xL))/(xU - xL));
Hadd(1,X.i) = -concv_slope;
Hadd(1,Z.i) = 1; 
kadd(1,1) = -concv_slope*xL + sec(xL);

% Convex lowerbound
% z > diff(sec(x),x)|(x=xL)*(x - xL) + sec(xL)
% z > diff(sec(x),x)|(x=xmid)*(x - xmid) + sec(xmid)
% z > diff(sec(x),x)|(x=xU)*(x - xU) + sec(xU)
Hadd(2,X.i) = diffsec(xL);
Hadd(3,X.i) = diffsec(xmid);
Hadd(4,X.i) = diffsec(xU);
Hadd(2,Z.i) = -1;
Hadd(3,Z.i) = -1;
Hadd(4,Z.i) = -1;
kadd(2,1) = diffsec(xL)*xL - sec(xL);
kadd(3,1) = diffsec(xmid)*xmid - sec(xmid);
kadd(4,1) = diffsec(xU)*xU - sec(xU);


% Adds the inequality constraints to polyrelax Hrep
Polyrelax.addHrepHrow(Hadd,kadd);        


% 20-05-2024: First version
    
end

function out = diffsec(x)
out = sin(x)/cos(x)^2;
end
