classdef Strip
% CLASS STRIP
% Defines the Strip object, properties, and methods
%
% CONSTRUCTOR
%       Z = Strip(p,d,sigma); - Creates a strip satisfying {|p'*x - d|<= sigma}
%       Z = Strip(p);         - Creates a strip satisfying {|p'*x|<= 1}
%       Z = Strip(p,d);       - Creates a strip satisfying {|p'*x - d|<= 1}
%
% PROPERTIES
%       p, d, sigma = defined according to the set definition
%
% METHODS
%       see help of each method

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
       
    properties (SetAccess=private)
        % Strip variables
        p;
        d;
        sigma;     
    end
     
    methods

        function obj = Strip(varargin)
            if nargin == 1 % {|p'*x|<= 1}
                p = varargin{1};
                if isnumeric(p)
                    if(size(p,2)==1)
                        obj.p = p;
                    else
                        error('Wrong input dimensions in Strip constructor.');
                    end
                else
                    error('Invalid input in Strip constructor.');
                end
                obj.d = 0;
                obj.sigma = 1;                
            elseif nargin == 2 % {|p'*x - d|<= 1}
                p = varargin{1};
                if isnumeric(p)
                    if(size(p,2)==1)
                        obj.p = p;
                    else
                        error('Wrong input dimensions in Strip constructor.');
                    end
                else
                    error('Invalid input in Strip constructor.');
                end
                d = varargin{2};
                if isnumeric(d)
                    if(numel(d)==1)
                        obj.d = d;
                    else
                        error('Wrong input dimensions in Strip constructor.');
                    end
                else
                    error('Invalid input in Strip constructor.');
                end
                obj.sigma = 1;
            elseif nargin == 3 % {|p'*x - d|<= sigma}
                p = varargin{1};
                if isnumeric(p)
                    if(size(p,2)==1)
                        obj.p = p;
                    else
                        error('Wrong input dimensions in Strip constructor.');
                    end
                else
                    error('Invalid input in Strip constructor.');
                end
                d = varargin{2};
                if isnumeric(d)
                    if(numel(d)==1)
                        obj.d = d;
                    else
                        error('Wrong input dimensions in Strip constructor.');
                    end
                else
                    error('Invalid input in Strip constructor.');
                end
                sigma = varargin{3};
                if ((isnumeric(sigma))&&(sigma>=0))
                    if(numel(sigma)==1)
                        obj.sigma = sigma;
                    else
                        error('Wrong input dimensions in Strip constructor.');
                    end
                else
                    error('Invalid input in Strip constructor.');
                end
            else
                error('The Strip constructor accepts only one, two or three input arguments')
            end
        end
        
        % Convenience function for dimension
        function val = dim(obj)
            val = size(obj.p,1);
        end
  
    end
    
end