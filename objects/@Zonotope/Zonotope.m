classdef Zonotope
% CLASS ZONOTOPE
% Defines the Zonotope object, properties, and methods
%
% CONSTRUCTOR
%       Z = Zonotope(c,G); - Creates a zonotope with center 'c' and generator matrix 'G'    
%       Z = Zonotope(a); - if 'a' is: - real vector: converts point 'a' into G-rep
%                                     - interval (INTLAB): converts the box 'a' into G-rep    
%       Z = Zonotope;      - Creates a zero-dimensional zonotope with no generators
%
% PROPERTIES
%       c: center
%       G: generator matrix
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
        % Center
        c;
        % Generators
        G;      
    end
   
    methods

        function obj = Zonotope(varargin)
            if nargin == 0 % Returns a null zonotope (dimension zero)
                obj = Zonotope.null;
            elseif nargin == 1 %Have only one entry
                c = varargin{1};
                if isa(c,'Zonotope') % If its already a zonotope, returns the input
                    obj = c;                                    
                elseif (isa(c,'intval')||isa(c,'Interval')) % It's an interval 
                    obj = Zonotope.frominterval(c);   
                elseif isnumeric(c) % its a point
                    if((size(c,2)==1))
                        dim = size(c,1);
                        obj.c = c;
                        obj.G = zeros(dim,0);
                    elseif(size(c,1)+size(c,2)==0) % [] input
                        obj = Zonotope.null;
                    else
                        error('Wrong input dimensions in Zonotope constructor.');
                    end
                else
                    error('Invalid input in Zonotope constructor.');
                end
            elseif nargin == 2 % Two inputs, c and G
                c = varargin{1};
                G = varargin{2};
                if(isnumeric(c)&&isnumeric(G))
                    if((size(c,2)==1)&&(size(c,1)==size(G,1)))
                        obj.c = c;
                        obj.G = G;
                    else
                        error('Wrong input dimensions in Zonotope constructor.');
                    end
                else
                    error('Invalid input in Zonotope constructor.');
                end
            else
                error('The Zonotope constructor accepts only zero, one or two input arguments')
            end
        end
        
        
        % Convenience functions for dimension and number of generators
        function val = dim(obj)
            val = size(obj.c,1);
        end
        function val = ng(obj)
            val = size(obj.G,2);
        end       
      
    end
    
    methods(Static)
        
        function Z = inclusion(c,M)
        %INCLUSION returns a zonotope containing a family of zonotopes (Theorem 3
        %          from Alamo et al, 2005)
        %
        %   SYNTAX: Z = INCLUSION(c,M)
        %
        %   INPUTS
        %           p: center common to the family of zonotopes
        %           M: interval matrix of generators (intval or interval object)
        %
        %   OUTPUT
        %           Z: zonotope containing the family of zonotopes

            diamM = diam(M);
            sumDiamM = sum(diamM,2);
            G = (1/2).*diag(sumDiamM);

            Z = Zonotope(c, [mid(M), G]);

        end
       
    end
    
    methods(Static, Access=private)
        
        function val = null
            val = Zonotope(zeros(0,1));
        end
        
        function val = frominterval(Z)
            val = Zonotope(mid(Z),diag(rad(Z)));
        end        
    end
end