classdef LZonotope
% CLASS LZONOTOPE
% Defines the line zonotope (LZonotope) object, properties, and methods
%
% CONSTRUCTOR
%       Z = LZonotope(c,G,M,S,A,b); - Creates a line zonotope with 'c', generators 'G', lines 'M', and equality constraint variables '(S,A,b)'
%       Z = LZonotope(c,[],M);      - Creates a line zonotope with center 'c' and lines 'M'     
%       Z = LZonotope(c,G,M);       - Creates a line zonotope with center 'c', generators 'G' and lines 'M'
%       Z = LZonotope(c,G,A,b);     - Creates a line zonotope with center 'c', generators 'G', and equality constraint variables '(A,b)'
%       Z = LZonotope(c,G);         - Creates a line zonotope with center 'c' and generators 'G'
%       Z = LZonotope(a); - if 'a' is: - real vector: converts 'a' into CLG-rep with center 'a'
%                                      - interval: converts the box 'a' into CLG-rep    
%                                      - zonotope: converts the zonotope 'a' into CLG-rep
%                                      - constrained zonotope: converts the constrained zonotope 'a' into CLG-rep    
%                                      - strip: converts the strip 'a' into CLG-rep       
%       Z = LZonotope;              - Creates a zero-dimensional line zonotope with no generators or lines
%
% PROPERTIES
%       c: center
%       G: generator matrix
%       M: lines matrix
%   S,A,b: equality constraints variables
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
        % Lines and generators
        G; M;
        % Equality constraints variables
        S; A; b;        
    end
   
    methods

        function obj = LZonotope(varargin)
            if nargin == 0 % Returns a null czonotope (dimension zero)
                obj = LZonotope.null;
            elseif nargin == 1 %Have only one entry
                c = varargin{1};
                if isa(c,'LZonotope') % If its already a LZonotope, returns the input
                    obj = c;                                    
                elseif (isa(c,'intval')||isa(c,'Interval')) % It's an interval 
                    obj = LZonotope.frominterval(c);   
                elseif isa(c,'Zonotope') % If its a zonotope, converts into CLG-rep
                    obj.c = c.c; 
                    obj.G = c.G;
                    obj.M = zeros(size(c.c,1),0);
                    obj.S = zeros(0,0);
                    obj.A = zeros(0,size(c.G,2));
                    obj.b = zeros(0,1);                    
                elseif isa(c,'CZonotope') % If its a constrained zonotope, converts into CLG-rep
                    obj.c = c.c; 
                    obj.G = c.G;
                    obj.M = zeros(size(c.c,1),0);
                    obj.S = zeros(size(c.A,1),0);
                    obj.A = c.A;
                    obj.b = c.b;
                elseif isa(c,'Strip') % If its a strip, converts into CLG-rep
                    strip_dim = size(c.p,1);
                    obj.c = zeros(strip_dim,1);
                    obj.G = zeros(strip_dim,1);
                    obj.M = eye(strip_dim);
                    obj.S = c.p.';
                    obj.A = -c.sigma;
                    obj.b = c.d;                    
                elseif isnumeric(c) % its a point
                    if((size(c,2)==1))
                        dim = size(c,1);
                        obj.c = c;
                        obj.G = zeros(dim,0);
                        obj.M = zeros(dim,0);
                        obj.A = zeros(0,0);                        
                        obj.S = zeros(0,0);
                        obj.b = zeros(0,1);                                                
                    elseif(size(c,1)+size(c,2)==0) % [] input
                        obj = LZonotope.null;
                    else
                        error('Wrong input dimensions in LZonotope constructor.');
                    end
                else
                    error('Invalid input in LZonotope constructor.');
                end
            elseif nargin == 2 % Two inputs, c and G
                c = varargin{1};
                G = varargin{2};
                if(isnumeric(c)&&isnumeric(G)) % Must be numeric variables
                    if((size(c,2)==1)&&(size(c,1)==size(G,1))) % Check input dimensions
                        obj.c = c;
                        obj.G = G;
                        obj.M = zeros(size(c,1),0);  
                        obj.S = zeros(0,0);
                        obj.A = zeros(0,size(G,2));
                        obj.b = zeros(0,1);                      
                    else
                        error('Wrong input dimensions in LZonotope constructor.');
                    end
                else
                    error('Invalid input in LZonotope constructor.');
                end
            elseif nargin == 3 % Three inputs, (c,G,M)
                c = varargin{1};
                G = varargin{2};
                M = varargin{3};
                if(isempty(G))
                    G = zeros(size(c,1),0);
                end
                if(isnumeric(c)&&isnumeric(G)&&isnumeric(M)) % Must be numeric variables
                    if((size(c,2)==1)&&(size(c,1)==size(G,1))&&(size(c,1)==size(M,1))) % Check input dimensions
                        obj.c = c;
                        obj.G = G;
                        obj.M = M;
                        obj.S = zeros(0,size(M,2));
                        obj.A = zeros(0,size(G,2));
                        obj.b = zeros(0,1);
                    else
                        error('Wrong input dimensions in LZonotope constructor.');
                    end
                else
                    error('Invalid input in LZonotope constructor.');
                end                
            elseif nargin == 4 % Four inputs, (c,G,A,b)
                c = varargin{1};
                G = varargin{2};
                A = varargin{3};
                b = varargin{4};
                if(isnumeric(c)&&isnumeric(G)&&isnumeric(A)&&isnumeric(b)) % Must be numeric variables
                    if((size(c,2)==1)&&(size(c,1)==size(G,1))&&(size(G,2)==size(A,2))&&(size(A,1)==size(b,1))&&(size(b,2)==1)) % Check input dimensions
                        obj.c = c;
                        obj.G = G;
                        obj.M = zeros(size(c,1),0);
                        obj.S = zeros(size(A,1),0);
                        obj.A = A;
                        obj.b = b;                      
                    else
                        error('Wrong input dimensions in LZonotope constructor.');
                    end
                else
                    error('Invalid input in LZonotope constructor.');
                end
            elseif nargin == 6 % Six inputs, (c,G,M,S,A,b)
                c = varargin{1};
                G = varargin{2};
                M = varargin{3};                
                S = varargin{4};                                
                A = varargin{5};
                b = varargin{6};
                if(isnumeric(c)&&isnumeric(G)&&isnumeric(M)&&isnumeric(S)&&isnumeric(A)&&isnumeric(b)) % Must be numeric variables
                    if((size(c,2)==1)&&(size(c,1)==size(G,1))&&(size(c,1)==size(M,1))...
                           &&(size(G,2)==size(A,2))&&(size(M,2)==size(S,2))&&(size(S,1)==size(b,1))&&(size(A,1)==size(b,1))&&(size(b,2)==1)) % Check input dimensions
                        obj.c = c;
                        obj.G = G;
                        obj.M = M;
                        obj.S = S;
                        obj.A = A;
                        obj.b = b;                      
                    else
                        error('Wrong input dimensions in LZonotope constructor.');
                    end
                else
                    error('Invalid input in LZonotope constructor.');
                end                
                
            else
                error('The LZonotope constructor accepts only zero, one, two, or four input arguments')
            end
        end        
        
        % Convenience functions for dimension and number of generators,
        % lines, and constraints
        function val = dim(obj)
            val = size(obj.c,1);
        end
        function val = ng(obj)
            val = size(obj.G,2);
        end
        function val = nl(obj)
            val = size(obj.M,2);
        end        
        function val = nc(obj)
            val = size(obj.b,1);
        end    
        
    end
    
    methods(Static)
        
        % Conversion from H-rep
        function Z = fromhrep(varargin)
        %LZONOTOPE.FROMHREP converts a convex polytope expressed in the half-space
        %                   representation {x : Hx <= k, Ax = b} into a
        %                   line zonotope (Theorem 1 in the CZ paper)
        %
        %   SYNTAX: Z = LZONOTOPE.FROMHREP(H,k)
        %           Z = LZONOTOPE.FROMHREP(H,k,A,b)
        %
        %   INPUTS
        %           H: matrix H from {x : Hx <= k, Ax = b}
        %           k: vector k from {x : Hx <= k, Ax = b}
        %           A: matrix A from {x : Hx <= k, Ax = b}
        %           b: vector b from {x : Hx <= k, Ax = b}
        %   OUTPUT
        %           Z: line zonotope as a structure
        
        % Just call czonotope.fromhrep
        Z = CZonotope.fromhrep(varargin);
        
        end
        
        function Z = realset(dim)
           % Convenience function returning the n-dim real set in CLG-rep
           Z = LZonotope(zeros(dim,1),zeros(dim,0),eye(dim));
        end
    end

    
    methods(Static, Access=private)
        
        function val = null
            val = LZonotope(zeros(0,1));
        end
        
        function val = frominterval(Z)
            val = LZonotope(mid(Z),diag(rad(Z)));
        end        
    end
end



