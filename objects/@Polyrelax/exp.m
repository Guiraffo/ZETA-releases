function Z = exp(X)
%EXP returns the exponential of a polyrelax object
%
%   SYNTAX: Z = EXP(X)
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

X = Polyrelax(X);
Z = Polyrelax(exp(X.x));


% Compute the envelopes

xL = inf(X.x);
xU = sup(X.x);
% xmid = 0.5*(xL+xU);
% Hadd = zeros(4,Z.i);
% kadd = zeros(4,1);

% Concave upperbound
% z < ((exp(xU) - exp(xL))/(xU - xL))*(x - xL) + exp(xL)
Hconcv = zeros(1,Z.i);
kconcv = zeros(1,1);
concv_slope = ((exp(xU) - exp(xL))/(xU - xL));
Hconcv(1,X.i) = -concv_slope;
Hconcv(1,Z.i) = 1; 
kconcv(1,1) = -concv_slope*xL + exp(xL);

% Convex lowerbound
% z > diff(exp(x),x)|(x=xL)*(x - xL) + exp(xL)
% z > diff(exp(x),x)|(x=xmid)*(x - xmid) + exp(xmid)
% z > diff(exp(x),x)|(x=xU)*(x - xU) + exp(xU)
switch toolsettings.polyrelax_approxmode
    case 0
        Hconvx = zeros(3,Z.i);
        kconvx = zeros(3,1);
        xmid = 0.5*(xL+xU);        
        Hconvx(1,X.i) = exp(xL);
        Hconvx(2,X.i) = exp(xmid);
        Hconvx(3,X.i) = exp(xU);
        Hconvx(1,Z.i) = -1;
        Hconvx(2,Z.i) = -1;
        Hconvx(3,Z.i) = -1;
        kconvx(1,1) = exp(xL)*xL - exp(xL);
        kconvx(2,1) = exp(xmid)*xmid - exp(xmid);
        kconvx(3,1) = exp(xU)*xU - exp(xU);
    case 2
        Hconvx = zeros(2,Z.i);
        kconvx = zeros(2,1);        
        xlow = (2/3)*xL + (1/3)*xU;
        xupp = (1/3)*xL + (2/3)*xU;
        Hconvx(1,X.i) = exp(xlow);
        Hconvx(2,X.i) = exp(xupp);
        Hconvx(1,Z.i) = -1;
        Hconvx(2,Z.i) = -1;
        kconvx(1,1) = exp(xlow)*xlow - exp(xlow);
        kconvx(2,1) = exp(xupp)*xupp - exp(xupp);           
    otherwise
        error('Invalid polyrelax approximation mode.')
end
    
% Adds the inequality constraints to polyrelax Hrep
% Polyrelax.addHrepHrow(Hadd,kadd);
Polyrelax.addHrepHrow([Hconcv;Hconvx],[kconcv;kconvx]);

% 21-05-2024: added approximation modes

end