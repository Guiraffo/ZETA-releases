function [Z,pivot_record,flag] = PartialSolve(Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    libCZon.PartialSolve(Z);                                   %
% Description:                                                            %
% Input:       Z      - constrained zonotope                              %
% Output:      Z      - constrained zonotope                              %
%              flag   - status flag. Zero if successful.                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    FuncID = 'libCZon_PartialSolve';
    flag=0;
    %---------------------------------------------------------%


    
    %---------------------------------------------------------%
    %Check inputs
    [n  ng] = size(Z{2});
    [nc ~ ] = size(Z{3});
    
    
    % MODIFICATION
%     if (nc==0); return; end;
    %Check inputs
    if (nc==0 || nc>ng)
        disp(['Error in ' FuncID]);
        nc
        pivot_record = 1:ng;
        keyboard;
        return;
    end;
    %---------------------------------------------------------%
    
    
    
    %---------------------------------------------------------%
    %
    G=Z{2};
    Ab=[Z{3} Z{4}];
    
    
    tol = max(nc,ng+1)*eps('double')*norm(Ab,'inf');
    
    % ORIGINAL CODE
%     pivot_record=1:ng;

    % CODE FROM PARTIALSOLVE
    col_pivots=1:ng+1;
    row_pivots=1:nc;
    for row=1:nc
    
        %Find the best pivot
        [~,xi_elim] = max(abs(Ab(row:nc,1:ng)),[],2);
        Ab_Scaled=Ab;
        Ab_Norms=ones(nc,1);
        
%      CODE MODIFIED ACCORDING TO PARTIALSOLVE.M, NOW INCLUDES TOLERANCE
%      CHECKING        
        for i=row:nc
            if ( abs(Ab(i,xi_elim(i-row+1)))>tol ) 
                Ab_Scaled(i,:) = Ab(i,:)/Ab(i,xi_elim(i-row+1));
                Ab_Norms(i) = norm(Ab_Scaled(i,:),1)-1;
            else
                %Everything in row i from 1:nc is zero
                Ab(i,1:nc)=0;
                Ab_Scaled(i,1:nc)=0;
                Ab_Norms(i)=Inf;  
            end
        end        
        
        % ORIGINAL CODE
%         for i=row:nc
%             Ab_Scaled(i,:) = Ab(i,:)/Ab(i,xi_elim(i-row+1));
%             Ab_Norms(i) = norm(Ab_Scaled(i,:),1)-1;
%         end
        
        
        [~,index] = sort(Ab_Norms(row:nc));
        pivot_row = index(1)+row-1;
        pivot_col = xi_elim(index(1));
        
        % ORIDINAL CODE
        %Permute
        %Ab([row pivot_row],:) = Ab([pivot_row row],:);
        %Ab(:, [row pivot_col]) = Ab(:, [pivot_col row]);
        %G(:, [row pivot_col]) = G(:, [pivot_col row]);
        %pivot_record([row pivot_col]) = pivot_record([pivot_col row]);
        
        % CODE FROM PARTIALSOLVE
        %Permute
        Ab([row pivot_row],:) = Ab([pivot_row row],:);
        row_pivots([row pivot_row]) = row_pivots([pivot_row row]);
        Ab(:, [row pivot_col]) = Ab(:, [pivot_col row]);
        G(:, [row pivot_col]) = G(:, [pivot_col row]);        
        col_pivots([row pivot_col]) = col_pivots([pivot_col row]);        
        
        
%      CODE MODIFIED ACCORDING TO PARTIALSOLVE.M (WITHOUT
%      RETURNING WHEN FLAG -4 )
        %Eliminate
        %if (Ab(row,row)~=0)
        if (abs(Ab(row,row))>tol)
            Ab(row,row:ng+1)=Ab(row,row:ng+1)/Ab(row,row);

            %Eliminate downward        
            for i=row+1:nc
                %if (abs(Ab(i,row))>1e-10)
                if (abs(Ab(i,row))>tol)
                    %Ab(i,:) = (1/Ab(i,row))*Ab(i,:)-Ab(row,:);
                    Ab(i,:) = Ab(i,:)-Ab(i,row)*Ab(row,:);
                end
                Ab(i,row)=0;
            end  

            %Eliminate upward
            for i=row-1:-1:1
                %if (abs(Ab(i,row))>1e-10)
                if (abs(Ab(i,row))>tol)
                    Ab(i,row:ng+1) = Ab(i,row:ng+1)-Ab(i,row)*Ab(row,row:ng+1);
                end
                Ab(i,row)=0;
            end   

        else
            % REMOVED THE RETURN STATEMENT
             disp('flag = -4 in libCZon_PartialSolve2');
        end
        
        
%      CODE USED IN THE RECENT AUTOMATICA PAPERS        
%         %Eliminate
%         %if (Ab(row,row)~=0)
%         if (abs(Ab(row,row))>tol)
%             Ab(row,row:ng+1)=Ab(row,row:ng+1)/Ab(row,row);
%         else
%            %flag=-4; return; 
%             disp('flag = -4 in libCZon_PartialSolve2');
%             %Z{2}=G;
%             %Z{3}=Ab(:,1:ng);
%             %Z{4}=Ab(:,ng+1);           
%            return
%         end
%         
%         %Eliminate downward        
%         for i=row+1:nc
%             %if (abs(Ab(i,row))>1e-10)
%             if (abs(Ab(i,row))>tol)
%                 Ab(i,:) = (1/Ab(i,row))*Ab(i,:)-Ab(row,:);
%             end
%             Ab(i,row)=0;
%         end  
%         
%         %Eliminate upward
%         for i=row-1:-1:1
%             %if (abs(Ab(i,row))>1e-10)
%             if (abs(Ab(i,row))>tol)
%                 Ab(i,row:ng+1) = Ab(i,row:ng+1)-Ab(i,row)*Ab(row,row:ng+1);
%             end
%             Ab(i,row)=0;
%         end
        
    end
    
    Z{2}=G;
    Z{3}=Ab(:,1:ng);
    Z{4}=Ab(:,ng+1);
    %---------------------------------------------------------%
    %pivot_record
    %pivot_record = col_pivots;
    %---------------------------------------------------------%
    pivot_record = col_pivots(1:ng);    % There is an extra column in these records
    
end