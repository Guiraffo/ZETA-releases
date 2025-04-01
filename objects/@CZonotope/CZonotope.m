classdef CZonotope
% CLASS CZONOTOPE
% Defines the constrained zonotope (CZonotope) object, properties, and methods
%
% CONSTRUCTOR
%       Z = CZonotope(c,G,A,b); - Creates a constrained zonotope with center 'c', generator matrix 'G', and equality constraint variables 'A' and 'b'    
%       Z = CZonotope(c,G);     - Creates a constrained zonotope with center 'c' and generator matrix 'G'
%       Z = CZonotope(a); - if 'a' is: - real vector: converts 'a' into CG-rep with center 'a'
%                                      - interval: converts the box 'a' into CG-rep    
%                                      - zonotope: converts the zonotope in G-rep into CG-rep
%       Z = CZonotope;          - Creates a zero-dimensional constrained zonotope with no generators
%
% PROPERTIES
%       c: center
%       G: generator matrix
%     A,b: equality constraints variables
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
        %Center
        c;
        %Generators
        G;
        %Equality constraints variables
        A; b;       
        
    end
   
    methods

        function obj = CZonotope(varargin)
            if nargin == 0 % Returns a null CZonotope (dimension zero)
                obj = CZonotope.null;
            elseif nargin == 1 %Have only one entry
                c = varargin{1};
                if isa(c,'CZonotope') % If its already a CZonotope, returns the input
                    obj = c;                                    
                elseif (isa(c,'intval')||isa(c,'Interval')) % It's an interval 
                    obj = CZonotope.frominterval(c);   
                elseif isa(c,'Zonotope') % If its a zonotope, converts into CG-rep
                    obj.c = c.c; 
                    obj.G = c.G;
                    obj.A = zeros(0,size(c.G,2));
                    obj.b = zeros(0,1);
                elseif isnumeric(c) % its a point
                    if((size(c,2)==1))
                        dim = size(c,1);
                        obj.c = c;
                        obj.G = zeros(dim,0);
                        obj.A = zeros(0,0);                        
                        obj.b = zeros(0,1);                                                
                    elseif(size(c,1)+size(c,2)==0) % [] input
                        obj = CZonotope.null;
                    else
                        error('Wrong input dimensions in CZonotope constructor.');
                    end
                else
                    error('Invalid input in CZonotope constructor.');
                end
            elseif nargin == 2 % Two inputs, c and G
                c = varargin{1};
                G = varargin{2};
                if(isnumeric(c)&&isnumeric(G)) % Must be numeric variables
                    if((size(c,2)==1)&&(size(c,1)==size(G,1))) % Check input dimensions
                        obj.c = c;
                        obj.G = G;
                        obj.A = zeros(0,size(G,2));
                        obj.b = zeros(0,1);                      
                    else
                        error('Wrong input dimensions in CZonotope constructor.');
                    end
                else
                    error('Invalid input in CZonotope constructor.');
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
                        obj.A = A;
                        obj.b = b;                      
                    else
                        error('Wrong input dimensions in CZonotope constructor.');
                    end
                else
                    error('Invalid input in CZonotope constructor.');
                end
                
            else
                error('The CZonotope constructor accepts only zero, one, two, or four input arguments')
            end
        end
        
        
        % Convenience functions for dimension and number of generators and
        % constraints
        function val = dim(obj)
            val = size(obj.c,1);
        end
        function val = ng(obj)
            val = size(obj.G,2);
        end
        function val = nc(obj)
            val = size(obj.A,1);
        end          
        
    end
    
    methods(Static)
        
        % Conversion from H-rep
        function Z = fromhrep(varargin)
        %CZONOTOPE.FROMHREP converts a convex polytope expressed in the half-space
        %                   representation {x : Hx <= k, Ax = b} into a constrained
        %                   zonotope (Theorem 1 in the CZ paper)
        %
        %   SYNTAX: Z = CZONOTOPE.FROMHREP(H,k)
        %           Z = CZONOTOPE.FROMHREP(H,k,A,b)
        %
        %   INPUTS
        %           H: matrix H from {x : Hx <= k, Ax = b}
        %           k: vector k from {x : Hx <= k, Ax = b}
        %           A: matrix A from {x : Hx <= k, Ax = b}
        %           b: vector b from {x : Hx <= k, Ax = b}
        %   OUTPUT
        %           Z: constrained zonotope as a structure

        % Input check
        if(nargin<2)
            error('Not enough input arguments.');
        elseif(nargin==2)
            H = varargin{1};
            k = varargin{2};            
            if(isnumeric(H)&&isnumeric(k))
                setdim = size(H,2);                
                A = zeros(0,setdim);
                b = zeros(0,1);        
                if((size(H,1)~=size(k,1))||(size(k,2)~=1))
                    error('Wrong input dimensions in CZonotope.fromhrep.')
                end
            else
                error('Invalid inputs in CZonotope.fromhrep.')
            end
        elseif(nargin==4)
            H = varargin{1};
            k = varargin{2};
            A = varargin{3};
            b = varargin{4};
            setdim = size(H,2);            
            if(~(isnumeric(H)&&isnumeric(k)&&isnumeric(A)&&isnumeric(b)))
                error('Invalid inputs in  CZonotope.fromhrep.');
            elseif((size(H,1)~=size(k,1))||(size(k,2)~=1)||(size(A,2)~=setdim)||(size(A,1)~=size(b,1))||(size(b,2)~=1))
                error('Wrong input dimensions in  CZonotope.fromhrep.')
            end
        else
            error('Invalid inputs in CZonotope.fromhrep.');
        end

        
        % Options
        OPTIONS.Display = 'off';

        space_dimension = size(H,2);
        
        % Z0 as the smallest box containing {x : Hx <= k, Ax = b} (Linear Programming)

        zetalower = zeros(space_dimension,1);
        zetaupper = zeros(space_dimension,1);

        for j=1:space_dimension

            f = zeros(space_dimension,1);
            f(j) = 1;
            Aineq = H;
            bineq = k;

            x = optim.solvelp( f,Aineq,bineq,A,b,[],[],OPTIONS); zetalower(j) =  x(j);
            x = optim.solvelp(-f,Aineq,bineq,A,b,[],[],OPTIONS); zetaupper(j) =  x(j);      

        end

        % Z0 in Grep
        c = (zetalower+zetaupper)/2;
        G = 0.5*diag(zetaupper-zetalower);

        % Lower bound of the interval [sigma,k] containing Hx, where x satisfies Hx <= k (sigma is given by the lower bound of the interval hull of H*Z0)
        sigma = H*c - sum(abs(H*G),2);

        % Intersection
        Zc = c;
        ZG = [G, zeros(size(G,1),size(H,1))];
        ZA = [H*G, diag(sigma-k)/2];
        Zb = (k+sigma)/2 - H*c;
        
        Z = CZonotope(Zc,ZG,ZA,Zb);

        end
    end

    
    methods(Static, Access=private)
        
        function val = null
            val = CZonotope(zeros(0,1));
        end
        
        function val = frominterval(Z)
            val = CZonotope(mid(Z),diag(rad(Z)),zeros(0,size(Z,1)),zeros(0,1));
        end        
    end
end