function Z = tan(X)
%TAN returns the tangent of a polyrelax object
%
%   SYNTAX: Z = TAN(X)
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

% Check if input is not contained in [3*pi/2,7*pi/2]
if((xL<=xmin)||(xL>=xmax)||(xU<=xmin)||(xU>=xmax))
    error(' The tangent function is currently implemented only for the interval [-pi/2,pi/2].');
end

Z = Polyrelax(tan(X.x)); 
    
if(xL >= 0) % Positive interval: convex function
        
    xmid = 0.5*(xL+xU);
    %xmid = atan(0.5*(tan(xL)+tan(xU)));
    Hadd = zeros(4,Z.i);
    kadd = zeros(4,1);

    % Concave upperbound
    % z < ((tan(xU) - tan(xL))/(xU - xL))*(x - xL) + tan(xL)
    concv_slope = ((tan(xU) - tan(xL))/(xU - xL));
    Hadd(1,X.i) = -concv_slope;
    Hadd(1,Z.i) = 1; 
    kadd(1,1) = -concv_slope*xL + tan(xL);

    % Convex lowerbound
    % z > diff(x^a,x)|(x=xL)*(x - xL) + xL^a
    % z > diff(x^a,x)|(x=xmid)*(x - xmid) + xmid^a
    % z > diff(x^a,x)|(x=xU)*(x - xU) + xU^a
    Hadd(2,X.i) = difftan(xL);
    Hadd(3,X.i) = difftan(xmid);
    Hadd(4,X.i) = difftan(xU);
    Hadd(2,Z.i) = -1;
    Hadd(3,Z.i) = -1;
    Hadd(4,Z.i) = -1;
    kadd(2,1) = difftan(xL)*xL - tan(xL);
    kadd(3,1) = difftan(xmid)*xmid - tan(xmid);
    kadd(4,1) = difftan(xU)*xU - tan(xU);

elseif(xU <= 0) % Negative interval: concave function

    xmid = 0.5*(xL+xU);
    Hadd = zeros(4,Z.i);
    kadd = zeros(4,1);

    % Concave upperbound
    % z < diff(x^a,x)|(x=xL)*(x - xL) + xL^a
    % z < diff(x^a,x)|(x=xmid)*(x - xmid) + xmid^a
    % z < diff(x^a,x)|(x=xU)*(x - xU) + xU^a
    Hadd(1,X.i) = -difftan(xL);
    Hadd(2,X.i) = -difftan(xmid);
    Hadd(3,X.i) = -difftan(xU);
    Hadd(1,Z.i) = 1;
    Hadd(2,Z.i) = 1;
    Hadd(3,Z.i) = 1;
    kadd(1,1) = -difftan(xL)*xL + tan(xL);
    kadd(2,1) = -difftan(xmid)*xmid + tan(xmid);
    kadd(3,1) = -difftan(xU)*xU + tan(xU);

    % Convex lowerbound        
    % z > ((tan(xU) - tan(xL))/(xU - xL))*(x - xL) + tan(xL)
    convx_slope = ((tan(xU) - tan(xL))/(xU - xL));
    Hadd(4,X.i) = convx_slope;
    Hadd(4,Z.i) = -1; 
    kadd(4,1) = convx_slope*xL - tan(xL);

else % Interval that contains zero

    % By finding the roots 
    % xCV and xCC by bisection method
    [xCV,xCC] = rootsbybisection(xL,xU);

    % Middle ("main") interval
    xLef = max(xL,xCC); % lower bound of main interval
    xRig = min(xU,xCV); % upper bound of main interval

    % Concave and convex secant slopes
    convx_slope = ((tan(xRig) - tan(xL))/(xRig - xL));
    concv_slope = ((tan(xU) - tan(xLef))/(xU - xLef));

    % Hrep of secants
    Hsec = zeros(2,Z.i); 
    ksec = zeros(2,1);

    % Convex secant
    % z > ((tan(xU) - tan(xL))/(xU - xL))*(x - xL) + tan(xL)
    Hsec(1,X.i) = convx_slope;
    Hsec(1,Z.i) = -1; 
    ksec(1,1) = convx_slope*xL - tan(xL);           

    % Concave secant
    % z < ((tan(xU) - tan(xL))/(xU - xL))*(x - xL) + tan(xL)
    Hsec(2,X.i) = -concv_slope;
    Hsec(2,Z.i) = 1; 
    ksec(2,1) = -concv_slope*xLef + tan(xLef);


    % Linear approximation for Right interval
    if(xU<=xCV) % If the Right interval does not exist
        HRig = zeros(0,Z.i);
        kRig = zeros(0,1);
    else % If it exists, linearize in 2 points (the first point would be equivalent to the convex secant)
        HRig = zeros(2,Z.i);
        kRig = zeros(2,1);               
        xRigM = 0.5*(xCV+xU);

        % Linearized convex lowerbound for Right interval
        % z > diff(x^a,x)|(x=xmid)*(x - xmid) + tan(xmid)
        % z > diff(x^a,x)|(x=xU)*(x - xU) + tan(xU)
        HRig(1,X.i) = difftan(xRigM);
        HRig(2,X.i) = difftan(xU);
        HRig(1,Z.i) = -1;
        HRig(2,Z.i) = -1;
        kRig(1,1) = difftan(xRigM)*xRigM - tan(xRigM);
        kRig(2,1) = difftan(xU)*xU - tan(xU);
    end

    % Linear approximation for Left interval
    if(xL>=xCC) % If the Left interval does not exist
        HLef = zeros(0,Z.i);
        kLef = zeros(0,1);
    else % If it exists, linearize in 2 points (the last, third point would be equivalent to the concave secant)
        HLef = zeros(2,Z.i);
        kLef = zeros(2,1);               
        xLefM = 0.5*(xL+xCC);

        % Linearized concave upper bound for Left interval
        % z < diff(x^a,x)|(x=xL)*(x - xL) + xL^a
        % z < diff(x^a,x)|(x=xmid)*(x - xmid) + xmid^a
        HLef(1,X.i) = -difftan(xL);
        HLef(2,X.i) = -difftan(xLefM);
        HLef(1,Z.i) = 1;
        HLef(2,Z.i) = 1;
        kLef(1,1) = -difftan(xL)*xL + tan(xL);
        kLef(2,1) = -difftan(xLefM)*xLefM + tan(xLefM);
    end    


    % Get the polyhedral enclosure Qj
    Hadd = [Hsec; HRig; HLef];
    kadd = [ksec; kRig; kLef];

end

% Adds the inequality constraints to polyrelax Hrep
Polyrelax.addHrepHrow(Hadd,kadd);        
    
end

function [xCV,xCC] = rootsbybisection(xL,xU)
% Find roots using fzero

persistent OPTIONS;
if(isempty(OPTIONS))
    %OPTIONS = optimset('Display','off');
    OPTIONS = optimset('fzero');
    OPTIONS.MyTOL = 1e-4;
end

bnds = max(abs([xL,xU]));

xCV = fzero(@(x) trigconstrCV(x,xL), [-OPTIONS.MyTOL,OPTIONS.MyTOL+bnds], OPTIONS);
xCC = fzero(@(x) trigconstrCC(x,xU), [-bnds-OPTIONS.MyTOL,OPTIONS.MyTOL], OPTIONS);

end
function out = difftan(x)
out = (sec(x))^2;
end
function out = trigconstrCV(x,xL)
% Implements secant = derivative for xCV
out = (tan(x)-tan(xL))/(x-xL) - difftan(x);
end
function out = trigconstrCC(x,xU)
% Implements secant = derivative for xCC
out = (tan(xU)-tan(x))/(xU-x) - difftan(x);
end