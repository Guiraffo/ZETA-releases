function Znew = propagate(Z,fname,OPTIONS)
%PROPAGATE propagates an interval through nonlinear function using interval
%          arithmetic         
%
%   SYNTAX: Znew = propagate(Z,fname)
%           Znew = propagate(Z,fname,OPTIONS)
%
%   INPUTS
%           B: interval object to be propagated
%       fname: function name of the nonlinear function
%     OPTIONS: (optional) OPTIONS structure claimed using util.propagopts
%              (default value is util.propagopts('IANat'))
%
%   OUTPUT
%         Znew: interval enclosing the propagated Z through f

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

if(nargin<2)
    error('Not enough input arguments.');
elseif(nargin==2)
    OPTIONS = util.propagopts('IANat');
end

method = OPTIONS.method;
f = util.fhandle(fname);

switch method
    case 'IANat'
        Znew = propagate_IANat(Z,f);
    case 'IAMV'
        Znew = propagate_IAMV(Z,f);
    otherwise
        error('Invalid propagation method in Interval/propagate.');
end

end

%% Interval arithmetic - Natural interval extension
function Znew = propagate_IANat(Z,f)
    Znew = f.eval(Z);
end

%% Interval arithmetic - Mean value extension
function Znew = propagate_IAMV(Z,f)
    Znew = f.eval(mid(Z)) + f.Jacob(Z)*(Z-mid(Z));
end
