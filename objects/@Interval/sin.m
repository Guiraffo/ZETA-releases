function Z = sin(X)
%SIN returns the sine of an interval
%
%   SYNTAX: Z = SIN(X)
%
%   INPUTS
%           X: the input interval
%
%   OUTPUT
%           Z: the sine of the input interval

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

if(numel(X)>1)
    error('This function is only implemented for scalars.');
end  

xmin2 = 7*pi/2;

xL = inf(X);
xU = sup(X);

% % Check if input is not contained in [3*pi/2,7*pi/2]
% if((xL<xmin1)||(xL>xmin2)||(xU<xmin1)||(xU>xmin2))
%     error(' The sine function is currently implemented only for the interval [3*pi/2,7*pi/2].');
% end

% Get lower and upper bounds
zL = lowerbound(xL,xU,xmin2);
zU = -lowerbound(-xU,-xL,xmin2);

Z = Interval(zL,zU);

end

% Lower bound
function zL = lowerbound(xL,xU,xmin2)

% Compute auxiliary variables
n1 = getn1(xL);

% Verify the existence of the three intervals in pag 120 of Scott's thesis
flag = (xU - 2*(n1 - 1)*pi >= xmin2);

if(flag)
    zL = -1;
else
    zL = min(sin(xL),sin(xU));
end

end

% Auxiliary functions
function n1 = getn1(xL)
% Obtain n1 and n2 using the n(x) function
n1 = floor(nx(xL));
end
function out = nx(x)
% Evaluates the n(x) function
out = x/(2*pi) + (1/4);
end