function [Z,ng,flag] = ReduceOrderChisci( Z , zo )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    libZon.ReduceOrderChisci                                   %
% Description: Reduces the order of a zonotope by recursively applying    %
%              Chisci's rule for reducing n+1 generators to a             %
%              parallelotope.                                             %
% Input:       Z      - initial zonotope                                  %
%              zo     - desired order                                     %
% Output:      Z      - reduced zonotope                                  %
%              ng     - number of generators in reduced zonotope          %
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
    flag = 0;
    FuncID = 'libZon_ReduceOrderChisciPlus';
    %---------------------------------------------------------%


    
    %---------------------------------------------------------%
    %Check inputs
    [n,ng] = size(Z{2});
    if (n <= 0 || ng <= 0) ;   flag=-1; return; end;
    if (zo<1);                 flag=-1; return; end;
    %---------------------------------------------------------%
    
    
    
    %---------------------------------------------------------%
    %Number of generators to eliminate
    N_Elim = ng-zo*n;
    if (N_Elim<=0); return; end;
    
    %Find an invertible subset of generators
    [Dual,flag,~,pivot_record]=util.PartialSolve([Z{2} eye(n)],n,ng);
    if flag==-4
       %There is no invertible set; resort to simple reduction
       disp(['Warning: ',FuncID,' resorting to libZon_ReduceOrderSimple'])
       [Z,ng,flag] = libZon.ReduceOrderSimple( Z , zo );
       return;
    else
        if (flag); return; end;
    end
    
    %At this stage, n generators are in ParTope, while the rest are in
    %Ordered_Gen. The idea is to remove one generator at a time form
    %Ordered_Gen, add it to ParTope to generate a matrix with n+1
    %generators called ParTopePlus, and apply Chisci's method for optimally
    %enclosing ParTopePlus with a parallelotope. This updates ParTope, and
    %the process is repeated until N_Elim generators have been eliminated.
    
    %The loop below eliminates one generator from Ordered_Gen in each pass
    %by adding it to ParTope and reducing back to a parallelotope.
    for i=1:N_Elim
        
        %Which generator in Ordered_Gen should be added to ParTope for
        %reduction? Can argue that t is a good generator if column of
        %ReductionDesirabilityMeasure is either very small, very large, or
        %nearly a scaled unit vector. Simple hueristic below considers only
        %the first two cases.
        %[ReductionDesirabilityMeasure, CondInv] = linsolve(ParTope,Other_Gen);
        R = abs(Dual(:,n+1:ng));
        RDM = prod(ones(size(R))+R) - 1 - sum(R);
        
        %Find max and min element of ReductionDesirabilityMeasure
        [MinRDM,col_index] = min(RDM);
        
        %Remove generator
        Dual(:,n+col_index)=[];
        ng = ng-1;
        
        %Reduce n+1 gnerators in ParTopePlus = [ParTope tnp1] to n using
        %Chisci's parallelotope rule.
        Dual(:,n+1:ng+n) = diag(1./(1+R(1:n,col_index)))*Dual(:,n+1:ng+n);
        
    end
    
    %Reform zonotope generator matrix
    Z{2} = Dual(:,ng+1:ng+n)\Dual(:,1:ng);
    %---------------------------------------------------------%
  
end