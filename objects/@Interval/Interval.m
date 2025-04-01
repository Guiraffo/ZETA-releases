classdef Interval
% CLASS Interval
%
% Interval object to be used for basic interval arithmetic operations
%
% Ways of creating an interval:
%       x = Interval       % Create an interval [0,0]
%       x = Interval(a)    % Create an interval [a,a]
%       x = Interval(a,b)  % Create an interval [a,b]
% 
% Type help Interval/(function) for specific help.
%
% PROPERTIES: x.LB (left endpoint)
%             x.UB (right endpoint)
%
% UNARY OPERATIONS: inf            (lower bound)
%                   sup            (upper bound)
%                   abs            (absolute value)
%                   real           (real number representation of a degenerate interval)
%                   mid            (midpoint)
%                   width
%
% BINARY OPERATIONS: +  (plus)
%                    -  (minus)
%                    *  (mtimes)
%                    /  (mrdivide)
%                    == (eq)
%                    ~= (ne)
%                    <  (lt)
%                    hull       (interval hull)
%                    intersec   (intersection)
%                    union

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
   
    properties (SetAccess=private, GetAccess=private)
        
        % Endpoints        
        LB = 0;        
        UB = 0;
        
    end    
   
    methods

        function obj = Interval(firstarg,varargin)
            if nargin == 0
                obj.LB = 0;
                obj.UB = 0;
            elseif nargin == 1 %Have only one entry
                if isa(firstarg,'Interval') % It's an interval
                    obj = firstarg;
                else
                    if numel(firstarg)==0 % It is empty
                        obj = Interval.empty;
                    elseif numel(firstarg)==1 % It's a real number
                        obj.LB = firstarg;
                        obj.UB = firstarg;                   
                    else
                        dim1 = size(firstarg,1);
                        dim2 = size(firstarg,2);
                        obj = Interval.zeros(dim1,dim2);
                        % For some reason, using loops is faster
                        for i=1:dim1
                            for j=1:dim2
                                obj(i,j).LB = firstarg(i,j);
                                obj(i,j).UB = firstarg(i,j);
                            end
                        end
%                         values = num2cell(firstarg);
%                         [obj.LB] = values{:};
%                         [obj.UB] = values{:};
                    end
                end
            elseif nargin == 2
                lower = firstarg;
                upper = varargin{1};
                if (numel(lower)==1 && numel(upper)==1)
                    if(lower<=upper)
                        obj.LB = lower;
                        obj.UB = upper;                
                    else
                        %obj.LB = varargin{1};
                        %obj.UB = firstarg;
                        %disp('The input arguments were switched to proceed.')
                        error('The lower bound must be lower or equal than the upper bound.')
                    end
                elseif(ndims(lower) <= 2 && ndims(upper) <= 2)
                    %error('The Interval class constructor currently accepts only scalars as inputs')
                    if(ndims(lower)==ndims(upper) && all(size(lower)==size(upper)))
                        if(all(lower<=upper))
                            dim1 = size(lower,1);
                            dim2 = size(lower,2);                            
                            obj = Interval.zeros(dim1,dim2);
                            for i=1:dim1
                                for j=1:dim2
                                    obj(i,j).LB = lower(i,j);
                                    obj(i,j).UB = upper(i,j);
                                end
                            end                                         
                        else
                            error('The lower bound must be lower or equal than the upper bound.')
                            %obj.LB = varargin{1};
                            %obj.UB = firstarg;
                            %disp('The input arguments were switched to proceed.')
                            %disp('The second argument (upper bound) must be greater or equal that the first argument (lower bound). Thus, the input arguments were switched to proceed.')
                            %error('The second argument (upper bound) must be greater or equal that the first argument (lower bound).')
                        end
                    else
                        error('Dimension mismatch in the Interval constructor.');
                    end
                else   
                    error('The Interval object supports only up to 2D matrices.')
                end
            else
                error('The Interval constructor accepts only zero, one or two input arguments')
            end
        end
        
        function out = inf(obj)
            out = reshape(vertcat(obj.LB),size(obj));
        end
        
        function out = sup(obj)
            out = reshape(vertcat(obj.UB),size(obj));
        end        
      
    end
    
    methods(Static)


        function val = zeros(varargin)
            if nargin == 0
                % For zeros with no arguments
                val = Interval;
            elseif nargin == 1
                % Use property default values
                val(varargin{:},varargin{:}) = Interval;
            else
                % Use property default values
                val(varargin{:}) = Interval;
            end
        end
 
         
        function val = ones(varargin)
            if (nargin == 0)
                % For ones with no arguments
                val = Interval(1);
            elseif any([varargin{:}] <= 0)
                % For ones with any dimension <= 0   
                val = Interval.empty(varargin{:});
            else
                % For ones(m,n,...)
                % Replicate [1,1]
                val = repmat(Interval(1),varargin{:});
            end
        end
        
        function val = unit(varargin)
            if (nargin == 0)
                % For unit interval with no arguments
                val = Interval(-1,1);
            elseif any([varargin{:}] <= 0)
                % For unit with any dimension <= 0   
                val = Interval.empty(varargin{:});
            else
                % For unit(m,n,...)
                % Replicate [-1,1]
                val = repmat(Interval(-1,1),varargin{:});
            end   
        end
        
        function val = rand(varargin)
            if (nargin == 0)
                % For zeros('Color')
                val = Interval(rand);
            elseif any([varargin{:}] <= 0)
                % For zeros with any dimension <= 0   
                val = Interval.empty(varargin{:});
            else
                % For zeros(m,n,...,'Color')
                % Use property default values
                val = Interval(rand(varargin{:}));
            end
        end        
        
    end
end