function [Z,xiL_Score,xiU_Score,flag,xim_array,xir_array] = Scale(Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    libCZon.Scale(Z);                                          %
% Description: Tightens the fundamental hypercube based on constraints    %
%              and rescales Z accordingly. Procedure does not change the  %
%              set, but makes dualization less conservative. For best     %
%              results, execute libCZon.PartialSolve prior to this        %
%              function.                                                  %
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
    FuncID = 'libCZon_Scale';
    flag=0; xiL_Score=[]; xiU_Score=[];
    xim_array = zeros(size(Z{2},2),1); % Added by Brenner
    xir_array = ones(size(Z{2},2),1); % Added by Brenner
    %---------------------------------------------------------%


    
    %---------------------------------------------------------%
    %Check inputs
    [n  ng] = size(Z{2});
    [nc ~ ] = size(Z{3});
    
    if (nc==0)
        disp(strcat('No constraints on input to ',FuncID));
        return;
    end
    %---------------------------------------------------------%
    
    
    
    %---------------------------------------------------------%
    %Scale
    tol = max(nc,ng)*eps('double')*norm(Z{3},'inf');
    xiL_Score=ones(1,ng)*Inf;
    xiU_Score=xiL_Score;
    for k=1:2
    for row=1:nc
    
        %Tighten xi_i\in [-1,1]
        aT = Z{3}(row,:);
        b  = Z{4}(row)  ;
        aT_Norm=norm(aT,1);
        for col=1:ng

            if ( abs(aT(col))>tol)
                    
                % Bounds obtained by solving the constraint for 'col'-th generator
                aT_NormM=aT_Norm-abs(aT(col));
                xiL_Update = (b/aT(col))-abs(aT_NormM/aT(col)); % \rho^L
                xiU_Update = (b/aT(col))+abs(aT_NormM/aT(col)); % \rho^U
                
                % Intersection of the result with [-1,+1]
                xi_L=max(-1,xiL_Update);
                xi_U=min( 1,xiU_Update);
                
                % New radius and midpoint
                xi_r = 0.5*(xi_U-xi_L);
                xi_m = 0.5*(xi_U+xi_L);
                
                if (xiL_Update>=-1)
                    xiL_Score(col) = 0;
                else
                    if (xi_r>tol)
                        xiL_Score(col) = min(xiL_Score(col),(abs(xiL_Update)-1)/xi_r);
                    end
                end
                    
                if (xiU_Update<=1)
                    xiU_Score(col) = 0;
                else
                    if (xi_r>tol)
                        xiU_Score(col) = min(xiU_Score(col),(abs(xiU_Update)-1)/xi_r);
                    end             
                end

                if ( abs(xi_r-1)>tol )

                    %Check feasibility here
                    if (xi_r<-tol)
                       disp(strcat('Empty set generated in ',FuncID))
                       Z{1}=[];
                       Z{2}=[];
                       Z{3}=[];
                       Z{4}=[];
                       return;
                    else
                       xi_r=max(0,xi_r); 
                    end

                    Z{1}=Z{1}+Z{2}(:,col)*xi_m;
                    Z{4}=Z{4}-Z{3}(:,col)*xi_m;
                    Z{2}(:,col)=Z{2}(:,col)*xi_r;
                    Z{3}(:,col)=Z{3}(:,col)*xi_r;

                    xim_array(col) = xi_m; % Added by Brenner (stores every xi_m)
                    xir_array(col) = xi_r; % Added by Brenner (stores every xi_r)
                    
                    % Note: if xi_m = 0 and xi_r = 1 (no rescaling), this
                    % code is not executed, and the initial values on
                    % xim_array and xir_array are kept (0 and 1)
                
                end
                
            end
            
        end
        
            [Max,maxi] = max(abs(Z{3}(row,:)));
        if (Max>tol)
            Z{4}(row  ) = Z{4}(row  )/Z{3}(row,maxi); 
            Z{3}(row,:) = Z{3}(row,:)/Z{3}(row,maxi); 
        end
        
    end
    end
    %---------------------------------------------------------%
    
end