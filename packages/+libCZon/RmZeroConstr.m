function [Z, nc, flag, ZeroInd] = RmZeroConstr(Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    libCZon.RmZeroConstr(Z)                                    %
% Description: Removes zero constraints.                                  %
% Input:       Z      - constrained zonotope                              %
% Output:      Z      - modified constrained zonotope                     %
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
    flag  = 0;
    %---------------------------------------------------------%


    
    %---------------------------------------------------------%
    %Check inputs
    [n  ng] = size( Z{2} );
    [nc ~ ] = size( Z{3} );
    %if (nc<1); return; end;
    if (nc<1); ZeroInd = []; return; end; % MODIFIED BY BRENNER
    %---------------------------------------------------------%
    
    
    
    %---------------------------------------------------------%
    %Remove trivial constraints
    tol = max(nc,ng)*eps('double')*norm(Z{3},'inf');
    RowMaxs = max(abs(Z{3}),[],2);
    ZeroInd = find(RowMaxs<=ones(nc,1)*tol);
    
    for i=1:length(ZeroInd)
        if abs(Z{4}(ZeroInd(i)))>tol
            %The set is empty
            flag=-4;
            Z{1}=[];
            Z{2}=[];
            Z{3}=[];
            Z{4}=[];
            return;
        end
    end
    Z{3}(ZeroInd,:)=[];
    Z{4}(ZeroInd,:)=[]; % MODIFIED BY BRENNER (avoids an empty column vector to be mistakely turned into a empty row vector by MATLAB)
    
    %Update nc
    nc = size(Z{3},1);
    %---------------------------------------------------------%
  
end