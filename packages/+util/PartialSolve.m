function [A,flag,row_pivots,col_pivots] = PartialSolve(A,NR,NC)
% PARTIALSOLVE implements Gauss-Jordan elimination with full pivoting
% Inputs:       A      - nr-by-nc matrix, nr<nc                         
%               NR     - index of the last row able to pivot              
%               NC     - index of the last column able to pivot           
% Outputs:      A      - nr-by-nc matrix, nr<nc                          
%               flag   - status flag. Zero if successful.  

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

    %---------------------------------------------------------%
    %Initialize
    FuncID = 'PartialSolve';
    flag=0;
    [nr nc] = size(A);
    col_pivots=1:nc;
    row_pivots=1:nr;
    %---------------------------------------------------------%


    
    %---------------------------------------------------------%
    %Check inputs
    if (nr==0 || nr>nc || NC>nc || NR>nr)
        disp(['Error in ' FuncID]);
        return;
    end;
    %---------------------------------------------------------%
    
    
    
    %---------------------------------------------------------%
    %
    tol = max(nr,nc)*eps('double')*norm(A,'inf');
    for row=1:NR
    
        %Find the best pivot
        [~,xi_elim] = max(abs(A(row:NR,row:NC)),[],2);
        xi_elim=row-1+xi_elim;
        A_Scaled=A;
        A_Norms=ones(NR+1-row,1);
        for i=row:NR
            if ( abs(A(i,xi_elim(i-row+1)))>tol ) 
                A_Scaled(i,:) = A(i,:)/A(i,xi_elim(i-row+1));
                A_Norms(i) = norm(A_Scaled(i,:),1)-1;
            else
                %Everything in row i from 1:NC is zero
                A(i,1:NC)=0;
                A_Scaled(i,1:NC)=0;
                A_Norms(i)=Inf;
            end
        end
        [~,index] = sort(A_Norms(row:NR));
        pivot_row = index(1)+row-1;
        pivot_col = xi_elim(index(1));
        
        %Permute
        A([row pivot_row],:) = A([pivot_row row],:);
        row_pivots([row pivot_row]) = row_pivots([pivot_row row]);
        A(:, [row pivot_col]) = A(:, [pivot_col row]);
        col_pivots([row pivot_col]) = col_pivots([pivot_col row]);
        
        %Eliminate
        if ( abs(A(row,row))>tol)
            A(row,row:nc)=A(row,row:nc)/A(row,row);
            
            %Eliminate downward
            for i=row+1:nr
                if (abs(A(i,row))>tol)
                    A(i,:) = A(i,:)-A(i,row)*A(row,:);
                else
                    A(i,row)=0;
                end
            end 
            
            %Eliminate upward
            for i=row-1:-1:1
                if (abs(A(i,row))>1e-10)
                    A(i,row:nc) = A(i,row:nc)-A(i,row)*A(row,row:nc);
                else
                    A(i,row)=0;
                end
            end
        end
        
    end
    %---------------------------------------------------------%
   
end