classdef Polyrelax
% CLASS POLYRELAX
% Defines the Polyhedral relaxation (POLYRELAX) object, properties, and methods
%
% CONSTRUCTOR
%      X = POLYRELAX(x), creates a Polyrelax object from x, where x is a
%          numeric variable, interval, or another Polyrelax object (for
%          the latter it will just copy the input to the output) 
%
% OBJECT PROPERTIES
%       x: interval object
%       i: Polyrelax object index
%
% CLASS PROPERTIES (implemented through persistent variables in static methods)
%    Hrep: relaxation half-space representation {Hz <= k, Az = b}
%       Z: augmented interval variable z of the factorized function
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
        % Interval
        x;
        % Index
        i;
        
    end  
   
    methods
        
        function obj = Polyrelax(varargin)
            if nargin == 0 % Not enough inputs
                error('Not enough inputs in Polyrelax.');
            elseif nargin == 1
                x = varargin{1};
                
                if(isa(x,'Polyrelax')||isempty(x))
                    % If its a Polyrelax object, or if its an empty input, copies the input object (maintaining the former's index )
                    obj = x;
                else
                    if(isnumeric(x))
                        % If its a constant, creates a new Polyrelax object with index -1
                        obj.x = x;
                        obj.i = -1;
                        return;
                    elseif(isa(x,toolsettings.IA_class))
                        % If its an interval, creates a new Polyrelax object from that interval
                        obj.x = x;
                    else
                        error('Invalid inputs in Polyrelax.');
                    end

                    
                    if((Polyrelax.Hrep_exists)&&(Polyrelax.Z_exists))
                        % If is not the first Polyrelax object, assigns a new index and adds an extra column of zeroes to the Hrep                        
                        thisHrep = Polyrelax.Hrep;
                        Hsize = size(thisHrep.H);
                        Asize = size(thisHrep.A);
                        obj.i = Hsize(2)+1;
                        Polyrelax.Hrep(struct('H',[thisHrep.H,zeros(Hsize(1),1)],'k',thisHrep.k,'A',[thisHrep.A,zeros(Asize(1),1)],'b',thisHrep.b));
                        
                        thisZ = Polyrelax.Z;
                        Polyrelax.Z([thisZ; obj.x]);
                        
                    elseif((~Polyrelax.Hrep_exists)&&(~Polyrelax.Z_exists))
                        % If is the first Polyrelax object, assigns the index 1 and initializes the Hrep                        
                        obj.i = 1;
                        Polyrelax.Hrep(struct('H',zeros(0,1),'k',zeros(0,1),'A',zeros(0,1),'b',zeros(0,1)));
                        Polyrelax.Z(obj.x);
                    else
                        error('Hrep and Z variables in Polyrelax are inconsistent.');
                    end
                        
                end
            else
                error('Too many inputs in Polyrelax.');
            end
        end              
      
    end
    
    methods(Static)
        
        % Clear static variables, setting them empty
        function clear
            Polyrelax.Hrep([]);
            Polyrelax.Z([]);
            Polyrelax.elimind([]);
            Polyrelax.anyspecialcon(0);
        end            
        
        % Hrep handler. Get (call without input), set (call with input)
        function val = Hrep(newHrep)
            persistent Hrep_;
            if nargin
                Hrep_ = newHrep;
            end
            val = Hrep_;  
        end  
        
        % Verifies if Hrep exists (it is empty if it doesnt)
        function val = Hrep_exists
            val = ~isempty(Polyrelax.Hrep);
        end          
        
        
        % Adds new inequalities to Hrep
        function addHrepHrow(Hadd,kadd)
            thisHrep = Polyrelax.Hrep;
            Polyrelax.Hrep(struct('H',[thisHrep.H; Hadd],'k',[thisHrep.k; kadd],'A',thisHrep.A,'b',thisHrep.b));
        end
        
        % Adds new equalities to Hrep
        function addHrepArow(Aadd,badd)
            thisHrep = Polyrelax.Hrep;
            Polyrelax.Hrep(struct('H',thisHrep.H,'k',thisHrep.k,'A',[thisHrep.A; Aadd],'b',[thisHrep.b; badd]));
        end
        
        
        % Z handler. Get (call without input), set (call with input)
        function val = Z(newZ)
            persistent Z_;
            if nargin
                Z_ = newZ;
            end
            val = Z_;  
        end        
        
        % Verifies if Z exists (it is empty if it doesnt)
        function val = Z_exists
            val = ~isempty(Polyrelax.Z);
        end          
        
        
        % Adds a new interval to Z
        function addZnew(Znew)
            thisZ = Polyrelax.Z;
            Polyrelax.Z([thisZ;Znew]);
        end
        
    
        
        % elimind handler. Get (call without input), set (call with input)
        function val = elimind(newelimind)
            persistent elimind_;
            if nargin
                elimind_ = newelimind;
            end
            val = elimind_;  
        end   
        
        % Verifies if elimind exists (it is empty if it doesnt)
        function val = elimind_exists
            val = ~isempty(Polyrelax.elimind);
        end   
        
        % Adds a new index to elimind
        function addelimind(newindex)
            thisElimind = Polyrelax.elimind;
            Polyrelax.elimind([thisElimind;newindex]);
        end
        
        % special constraints flag handler. Get (call without input), set (call with input)
        function val = anyspecialcon(newvalue)
            persistent anyspecialcon_;
            if nargin
                anyspecialcon_ = newvalue;
            end
            val = anyspecialcon_;  
        end          
        
        % Solves all the equality constraints while transfering the
        % constraints defining variables of interest to the projector
        % matrix Eh and bias ch
        function [Gh,ch] = solveAbTriang(index_h)
            thisElimind = Polyrelax.elimind;
            thisHrep = Polyrelax.Hrep;
            H = thisHrep.H;
            k = thisHrep.k;
            A = thisHrep.A;
            b = thisHrep.b;
            nof_factors = size(thisHrep.H,2);
            nonelimind = setdiff(1:size(thisHrep.A,2),thisElimind);
                       
            % 'e' stands for factors to eliminate, 'n' for factors to not eliminate
            % He*ze + Hn*zn <= k;
            % Ae*ze + An*zn = b;
            He = H(:,thisElimind); Hn = H(:,nonelimind);            
            Ae = A(:,thisElimind); An = A(:,nonelimind);
            
            if(det(Ae)==0)
                error('Ae is not invertible in Polyrelax.solveAbTriang');
            end
            newH = Hn - (He/Ae)*An;
            newk = k - (He/Ae)*b;
            newA = zeros(0,length(nonelimind));
            newb = zeros(0,1);
            
            thisZ = Polyrelax.Z;
            newZ = thisZ(nonelimind,:);
            
            % Generate Gh and ch
            % First build the original projection matrix
            h_dim = length(index_h);
            Eh = zeros(h_dim,nof_factors);
            for j=1:h_dim
                Eh(j,index_h(j)) = 1;
            end
            % Reorder Eh as Ehe and Ehn
            Ehe = Eh(:,thisElimind);
            Ehn = Eh(:,nonelimind);
            % Build Gh and ch
            Gh = Ehn - (Ehe/Ae)*An;
            ch = (Ehe/Ae)*b;
            
            
            % Updates Hrep, Z, and elimind
            Polyrelax.Hrep(struct('H',newH,'k',newk,'A',newA,'b',newb));            
            Polyrelax.Z(newZ);
            Polyrelax.elimind([]);
        end   
        
        % Solves all the equality constraints resulting from linear
        % operations defining new factors, while transfering any
        % constraints defining the final factors to matrix Eh and bias ch.
        function [Gh,ch,inputind] = solveAbTriangPartial(index_h,nof_inps)
            inputind = 1:nof_inps;
            if(~Polyrelax.anyspecialcon) % If no special constraints are present, just call the old method
                [Gh,ch] = Polyrelax.solveAbTriang(index_h);
                return;
            end
            
            % --- Elimination of ze (triangular matrix) ------
            
            thisElimind = Polyrelax.elimind;
            thisHrep = Polyrelax.Hrep;
            H = thisHrep.H;
            k = thisHrep.k;
            A = thisHrep.A;
            b = thisHrep.b;
            nof_factors = size(thisHrep.H,2);
            nonelimind = setdiff(1:size(thisHrep.A,2),thisElimind);
            nof_elimind = length(thisElimind);
                       
            % 'e' stands for factors to eliminate, 'n' for factors to not eliminate
            % He*ze + Hn*zn <= k;
            % Ae*ze + An*zn = b;
            He = H(:,thisElimind); Hn = H(:,nonelimind);            
            Ae = A(:,thisElimind); An = A(:,nonelimind);
            
            % Find the rows that defined ze
            elimrows = zeros(1,nof_elimind);
            for j=1:nof_elimind
                elimrows(j) = find(A(:,thisElimind(j)),1); % Get the row associated to the first non-zero entry (this row defined ze)
            end
            non_elimrows = setdiff(1:size(thisHrep.A,1),elimrows);
            
            Aeu = Ae(elimrows,:);       Anu = An(elimrows,:);       bu = b(elimrows,:);
            Aed = Ae(non_elimrows,:);   And = An(non_elimrows,:);   bd = b(non_elimrows,:);
            
            if(det(Aeu)==0)
                error('Aeu is not invertible in Polyrelax.solveAb');
            end
            finalA = And - Aed*(Aeu\Anu);
            finalb = bd - Aed*(Aeu\bu);
            finalH = Hn - He*(Aeu\Anu);
            finalk = k - He*(Aeu\bu);            
            
            % Generate tilGh and tilch
            % First build the original projection matrix
            h_dim = length(index_h);
            Eh = zeros(h_dim,nof_factors);
            for j=1:h_dim
                Eh(j,index_h(j)) = 1;
            end
            % Reorder Eh as Ehe and Ehn
            Ehe = Eh(:,thisElimind);
            Ehn = Eh(:,nonelimind);
            % Build tilGh and tilch
            finalGh = Ehn - Ehe*(Aeu\Anu);
            finalch = Ehe*(Aeu\bu);
            
            % ------------------------------------------------            
            
            % Get Zf
            thisZ = Polyrelax.Z;
            finalZ = thisZ(nonelimind,:);            
            
            % Updates Hrep, Z, and elimind
            Polyrelax.Hrep(struct('H',finalH,'k',finalk,'A',finalA,'b',finalb));            
            Polyrelax.Z(finalZ);
            Polyrelax.elimind([]);
            
            % Output variables
            Gh = finalGh;
            ch = finalch;
        end                    
        
        
    end
end