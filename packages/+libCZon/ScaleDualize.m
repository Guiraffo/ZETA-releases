function [Z, ng, nc, flag, xi_elim, xi_m, xi_r, Gbar, Abar, bbar, pivot_row] = ScaleDualize(Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    libCZon.ScaleDualize(Z);                                   %
% Description:                                                            %
% Input:       Z      - constrained zonotope                              %
% Output:      Z     - dualized zonotope                                  %
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
    FuncID = 'libCZon_ScaleDualize';
    flag  = 0;
    xi_elim=[];
    %---------------------------------------------------------%


    
    %---------------------------------------------------------%
    %Get sizes
    [n  , ng] = size( Z{2} );
    [nc ~ ] = size( Z{3} );
    %---------------------------------------------------------%
    
    
    
    %---------------------------------------------------------%
    %Partial solve and scaling
%     [AG,flag] = PartialSolve([Z{3} -Z{4}; Z{2} Z{1}],nc,ng);
%     if (flag); 
%         disp(['PartialSolve failed in ',FuncID]); return; 
%     end;
%     Z{1} =  AG(nc+1:nc+n,ng+1);
%     Z{2} =  AG(nc+1:nc+n,1:ng);
%     Z{3} =  AG(1   :nc  ,1:ng);
%     Z{4} = -AG(1   :nc  ,ng+1);

%     %Partial solve and scaling
    [Z,~,flag] = libCZon.PartialSolve(Z);
    if (flag); 
        disp(['PartialSolve failed in ',FuncID]); return; 
    end;
    Gbar = Z{2};
    Abar = Z{3};
    bbar = Z{4};

    
    %Scale
    %[Z,xiL_Score,xiU_Score,flag] = libCZon_Scale(Z);
    [Z,xiL_Score,xiU_Score,flag,xi_m,xi_r] = libCZon.Scale(Z);       
    %[~,xiL_Score,xiU_Score,flag,~,~] = libCZon_Scale(Z);
    %[Z,xi_m,xi_r] = czonotope_rescaling(libCZon_Wrapfrom(Z),'LP');
    %Z = libCZon_Wrapto(Z);
    if (flag); disp(['libCZon_Scale failed in ',FuncID]); return; end;  
    
    %Remove trivial constraints
    nc0=nc;
    %[Z, nc, flag] = libCZon_RmZeroConstr(Z);
    [Z, nc, flag, ZeroInd] = libCZon.RmZeroConstr(Z); % MODIFIED BY BRENNER
    if (flag)
        if (flag==-4)
            disp(['Empty set identified in ',FuncID]); return;
        else
            disp(['libCZon_RmZeroConstr failed in ',FuncID]); return;
        end
    end
    
    % ADDED BY BRENNER IN 06-02-2019
    OriginalInd = 1:nc0;    
    if(nc-nc0 < 0)
       %disp(['Null constraints were removed in ',FuncID]);
       %pivot_row = ZeroInd;
       flag = -9;
       %return;
       OriginalInd(ZeroInd) = []; % Suppress the indexes of removed null constraints 
    else
       ZeroInd = [];
    end
        
    %Number of constraints left to dualize
    N_Dual = 1-(nc0-nc);
    N_Dual = min(N_Dual,nc);
    %if (N_Dual<=0); return; end;
    if (N_Dual<=0); flag = -9; pivot_row = [ZeroInd; -1]; return; end; % CODE MODIFIED BY BRENNER in 06-02-2019
    %---------------------------------------------------------%
    
    
    
    %---------------------------------------------------------%
    %Find best xi to eliminate
    
    xi_FScore = max(xiL_Score,xiU_Score);
    [FScore,F_xi_elim]=min(xi_FScore);
    if (FScore<1e-3)
        xi_elim = F_xi_elim;
    else
        
        %Need to compute G_Scores
        G = Z{2};
        A=Z{3};
        nc=size(A,1);
        ei=zeros(ng+nc,1);
        %Atil=[A; ei'];
        
        
        MMsub   = [(G'*G)+eye(ng) A'; A zeros(nc)];
        Mbsub   = [zeros(ng+nc,1)];
        [L,U,p] = lu(MMsub,'vector'); 
        Mx=zeros(ng,1);
        
        
        %Used for the old way below
        %MM = [MMsub zeros(ng+nc,1); zeros(1,ng+nc) 0];
        %Mb = [Mbsub; 0];
        
        
        for i=1:ng
            
            if (xi_FScore(i)==Inf)
                xi_GScore(i)=Inf;
            else
                
                ei(i)=1;
                
                opts.LT=true;
                uhat = linsolve(L,ei(p),opts);
                opts.LT=false;
                
                opts.UT=true;
                Ihat = linsolve(U,uhat,opts);
                opts.UT=false;
                
                Mx=Ihat(1:ng)*(xi_FScore(i)/Ihat(i));
                
                xi_GScore(i) = norm(G*Mx,2);
                
                ei(i)=0;
                
                %This is the old way to calculate the same thing
                %MM(ng+nc+1,i)=1;
                %MM(i,ng+nc+1)=1;
                %Mb(ng+nc+1)=xi_FScore(i);
                
                %opts.SYM=true;
                %[Mx1, CondRecip] = linsolve(MM,Mb,opts);
                %opts.SYM=false;
                %if (CondRecip<1e-12)
                %    disp('wft in '); disp(FuncID);
                %    flag=-1; return;
                %end;
                
                %error = Mx1(1:ng)-Mx
                
                %MM(ng+nc+1,i)=0;
                %MM(i,ng+nc+1)=0;
                
            end
            
        end
        [~,xi_elim]=min(xi_GScore);
        %xi_FScore=xi_FScore
        %xi_GScore=xi_GScore
    end
    %---------------------------------------------------------%
    
    
    
    %---------------------------------------------------------%
    %Eliminate
%    if (xi_elim>nc)
        [max_val,pivot_row]=max(abs(Z{3}(:,xi_elim)));
        if (max_val==0)
            disp(['Zero pivot in ',FuncID]);
            flag=-1; return;
        end
        
        
        % Complexity: O(nc(1+ng))
        
        Z{4}(pivot_row)   = Z{4}(pivot_row  )/Z{3}(pivot_row,xi_elim);
        Z{3}(pivot_row,:) = Z{3}(pivot_row,:)/Z{3}(pivot_row,xi_elim);
        
        for row=1:pivot_row-1 
            lambda = Z{3}(row,xi_elim);
            Z{4}(row  ) = Z{4}(row  )-lambda*Z{4}(pivot_row  ); 
            Z{3}(row,:) = Z{3}(row,:)-lambda*Z{3}(pivot_row,:); 
        end
        
        for row=pivot_row+1:nc
            lambda = Z{3}(row,xi_elim);
            Z{4}(row  ) = Z{4}(row  )-lambda*Z{4}(pivot_row  ); 
            Z{3}(row,:) = Z{3}(row,:)-lambda*Z{3}(pivot_row,:); 
        end
            
        for row=1:n % Complexity: O(n(1+ng))
            lambda = Z{2}(row,xi_elim);
            Z{1}(row  ) = Z{1}(row  )+lambda*Z{4}(pivot_row  ); 
            Z{2}(row,:) = Z{2}(row,:)-lambda*Z{3}(pivot_row,:); 
        end
        
%    else
%        pivot_row=xi_elim;
%    end
    
    Z{2}(:,xi_elim)=[];
    Z{3}(:,xi_elim)=[];
    Z{3}(pivot_row,:)=[];
    Z{4}(pivot_row,:)=[]; % MODIFIED BY BRENNER (avoids an empty column vector to be mistakely turned into a empty row vector by MATLAB)

    %Update ng and nc
    [nc,ng]=size(Z{3});
    
    
    % ADDED BY BRENNER IN 06-02-2019
    pivot_row = [ZeroInd; OriginalInd(pivot_row)];
    
    %---------------------------------------------------------%
    
    
end