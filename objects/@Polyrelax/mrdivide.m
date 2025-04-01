function Z = mrdivide(X,Y)
%/ or MRDIVIDE returns the product of two Polyrelax objects
%
%   SYNTAX: Z = X/Y
%           Z = MRDIVIDE(X,Y)
%
%   INPUTS
%           Z: the first Polyrelax object
%           W: the second Polyrelax object
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

X = Polyrelax(X);
Y = Polyrelax(Y);



% Adds four inequality constraints to Polyrelax Hrep associated to X = Y*Z
if(X.i==-1) % If X is a constant
    
    Z = Polyrelax(X.x/Y.x);
    
    x = X.x;
    yL = inf(Y.x);
    yU = sup(Y.x);
    zL = inf(Z.x);
    zU = sup(Z.x);

    Hadd = zeros(4,Z.i);
    Hadd(:,Y.i) = [zL; zU; -zU; -zL];
    Hadd(:,Z.i) = [yL; yU; -yL; -yU];
    kadd = [yL*zL + x; yU*zU + x; -yL*zU - x; -yU*zL - x];
    
elseif(Y.i==-1) % If Y is a constant
    
    Z = (1/Y.x)*X;
    return;
    
else % If neither X or Y is a constant, do x = y*z
    
    Z = Polyrelax(X.x/Y.x);

    yL = inf(Y.x);
    yU = sup(Y.x);
    zL = inf(Z.x);
    zU = sup(Z.x);

    Hadd = zeros(4,Z.i);
    Hadd(:,Y.i) = [zL; zU; -zU; -zL];
    Hadd(:,Z.i) = [yL; yU; -yL; -yU];
    Hadd(:,X.i) = [-1; -1;   1;   1];
    kadd = [yL*zL; yU*zU; -yL*zU; -yU*zL];

end


Polyrelax.addHrepHrow(Hadd,kadd);

% 20-04-2024: properly uses multiplication by a constant when y is a constant
% 05-03-2024: added proper treatment for constants
% 05-02-2024: first version
    
end