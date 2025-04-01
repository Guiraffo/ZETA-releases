function Z = mpower(X,a)
%^ or MPOWER returns the (positive integer) power of a interval
%
%   SYNTAX: Z = X^a
%           Z = MPOWER(X,a)
%
%   INPUTS
%           Z: the Interval object
%           a: the exponent (a positive integer)
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
if(~isnumeric(a))
    error('The second argument must be a numeric variable.');
end
if(numel(a)>1)
    error('The second argument must be a scalar.');
end
if(a <= 0)
    error('Only positive exponents are currently implemented.');
end

X = Interval(X);

% Call the respective subroutine
if(rem(a, 2)==1) % Odd exponent
    
Z = power_odd(X,a);
    
elseif(rem(a, 2)==0) % Even exponent

Z = power_even(X,a);    
    
else
%    error('Only integer exponents are implemented.');
Z = power_noninteger(X,a);    

end
    
end

function result = power_odd(X,a) % For odd integer power
    result = Interval(X.LB^a,X.UB^a);  
end

function result = power_even(X,a) % For even integer power
    endpts = [abs(X.LB),abs(X.UB)];
    if(isinside(X,0))
        result = Interval(0,max(endpts)^a);
    else
        result = Interval(min(endpts)^a,max(endpts)^a);  
    end
end

function result = power_noninteger(X,a) % For noninteger power (ONLY FOR X>=0)
    xL = inf(X);
    xU = sup(X);
    if(xL<=0)
        error('For non-integer power, only positive arguments are currently implemented.');
    end
    result = Interval(xL^a,xU^a);
end