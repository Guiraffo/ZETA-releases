function eq(X,a)
%== or EQUAL adds equality constraints between a column vector of polyrelax
%            objects and a real vector
%
%   SYNTAX: X==a
%           EQ(X,a)
%
%   INPUTS
%           X: the column vector of polyrelax objects
%           a: the real vector
%
%   OUTPUT
%           None

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

if(~isa(X,'Polyrelax')||~isnumeric(a))
    error('The first and second inputs must be a polyrelax object and a real variable, respectively.')
end

if((size(X,1)~=size(a,1))||(size(X,2)~=1)||(size(a,2)~=1))
   error('Wrong dimensions in the inputs.'); 
end

Hrep = Polyrelax.Hrep;
lastindex = size(Hrep.H,2);
Xinds = [X.i];
Xinds = Xinds(:); % Makes column vector for addelimind
nof_equalities = length(Xinds);
Aadd = zeros(nof_equalities,lastindex);

% Adds equality constraints to polyrelax Hrep
% [ ... 1 ... ]z = a
for j=1:nof_equalities
    Aadd(j,Xinds(j)) = 1;
end
badd = a;
   
Polyrelax.addHrepArow(Aadd,badd);
Polyrelax.anyspecialcon(1);

% 15-05-2024: Now sets the special constraint flag to 1
% 03-05-2024: Properly dealing with vectors instead of scalars
% 14-04-2024: Check if inputs are column vectors with same dimensions
% 12-04-2024: first version
    
end