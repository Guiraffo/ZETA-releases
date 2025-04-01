function Z = mtimes(X,Y)
%* or MTIMES returns the product of two Polyrelax objects
%
%   SYNTAX: Z = X*Y
%           Z = MTIMES(X,Y)
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

Z = Polyrelax(X.x*Y.x);


if(X.i==-1) % If X is a constant
    
    Aadd = zeros(1,Z.i);
    badd = zeros(1,1);    
    
    % [ ... -x ... 1]z = 0
    Aadd(Y.i) = -X.x;
    Aadd(Z.i) = 1;    
    badd = 0;
    
    Polyrelax.addHrepArow(Aadd,badd);
    Polyrelax.addelimind(Z.i);
    
elseif(Y.i==-1) % If Y is a constant
    
    Aadd = zeros(1,Z.i);
    badd = zeros(1,1);    
    
    % [ ... -x ... 1]z = 0
    Aadd(X.i) = -Y.x;
    Aadd(Z.i) = 1;    
    badd = 0;
    
    Polyrelax.addHrepArow(Aadd,badd);
    Polyrelax.addelimind(Z.i);
    
else % If neither X or Y is a constant
    % Adds four inequality constraints to Polyrelax Hrep from bilinearity
    % - McCormick relaxation

    xL = inf(X.x);
    xU = sup(X.x);
    yL = inf(Y.x);
    yU = sup(Y.x);
    %zL = inf(Z.x);
    %zU = sup(Z.x);

    Hadd = zeros(4,Z.i);
    Hadd(:,X.i) = [yL; yU; -yU; -yL];
    Hadd(:,Y.i) = [xL; xU; -xL; -xU];
    Hadd(:,Z.i) = [-1; -1;   1;   1];
    kadd = [xL*yL; xU*yU; -xL*yU; -xU*yL];
    
    Polyrelax.addHrepHrow(Hadd,kadd);    

end

% Revision 06-03-2024: multiplication by a constant leads to one equality constraint, instead of 4 inequalities
% Revision 05-03-2024: added proper treatment for constants
% First version: 05-02-2024
    
end