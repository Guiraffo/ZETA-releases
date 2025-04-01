function [Z, ng, flag] = RmZeroGen(Z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    libCZon.RmZeroGen(Z);                                      %
% Description: Removes zero (lifted) generators.                          %
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
    %---------------------------------------------------------%
    
    
    
    %---------------------------------------------------------%
    %Remove trivial generators
    ZeroGen = all(~[Z{2}; Z{3}],1);
    ZeroInd=find(ZeroGen);
    Z{2}(:,ZeroInd)=[];
    Z{3}(:,ZeroInd)=[];
    
    %Update ng
    ng = size(Z{2},2);
    %---------------------------------------------------------%
  
end