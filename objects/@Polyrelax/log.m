function Z = log(X)
%EXP returns the natural logarithm of a polyrelax object
%
%   SYNTAX: Z = LOG(X)
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
Z = Polyrelax(log(X.x));

if(~(inf(X.x)>0))
    error('Only positive scalars are allowed.')
end


% Compute the envelopes

xL = inf(X.x);
xU = sup(X.x);
% xmid = 0.5*(xL+xU);
% Hadd = zeros(4,Z.i);
% kadd = zeros(4,1);

% Concave upperbound
% z < diff(log(x),x)|(x=xL)*(x - xL) + log(xL)
% z < diff(log(x),x)|(x=xmid)*(x - xmid) + log(xmid)
% z < diff(log(x),x)|(x=xU)*(x - xU) + log(xU)
switch toolsettings.polyrelax_approxmode
    case 0
        Hconcv = zeros(3,Z.i);
        kconcv = zeros(3,1);
        xmid = 0.5*(xL+xU);        
        Hconcv(1,X.i) = -1/xL;
        Hconcv(2,X.i) = -1/xmid;
        Hconcv(3,X.i) = -1/xU;
        Hconcv(1,Z.i) = 1;
        Hconcv(2,Z.i) = 1;
        Hconcv(3,Z.i) = 1;
        kconcv(1,1) = -(1/xL)*xL + log(xL);
        kconcv(2,1) = -(1/xmid)*xmid + log(xmid);
        kconcv(3,1) = -(1/xU)*xU + log(xU);
    case 2
        Hconcv = zeros(2,Z.i);
        kconcv = zeros(2,1);
        xlow = (2/3)*xL + (1/3)*xU;
        xupp = (1/3)*xL + (2/3)*xU;
        Hconcv(1,X.i) = -1/xlow;
        Hconcv(2,X.i) = -1/xupp;
        Hconcv(1,Z.i) = 1;
        Hconcv(2,Z.i) = 1;
        kconcv(1,1) = -(1/xlow)*xlow + log(xlow);
        kconcv(2,1) = -(1/xupp)*xupp + log(xupp);        
    otherwise
        error('Invalid polyrelax approximation mode.')
end        
                

% Convex lowerbound
% z > ((log(xU) - log(xL))/(xU - xL))*(x - xL) + log(xL)
Hconvx = zeros(1,Z.i);
kconvx = zeros(1,1);
convx_slope = ((log(xU) - log(xL))/(xU - xL));
Hconvx(1,X.i) = convx_slope;
Hconvx(1,Z.i) = -1; 
kconvx(1,1) = convx_slope*xL - log(xL);

    
% Adds the inequality constraints to polyrelax Hrep
Polyrelax.addHrepHrow([Hconcv;Hconvx],[kconcv;kconvx]);

% 21-05-2024: added approximation modes

end