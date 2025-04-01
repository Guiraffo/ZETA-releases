function [Phi,Psi,Tildevars] = LTI_reachable_sets(r,A,Bu,Bw,s,C,Dv,X0,W,V,u)
%LTI_REACHABLE_SETS returns the reachable state and output sets of a linear
%                   time-invariant system given an initial set, disturbance
%                   bounds, and an input sequence (sets described as
%                   zonotopes)
%
%   SYNTAX: [Phi,Psi,Tildevars] = LTI_REACHABLE_SETS(A,B,Dw,r,C,Dv,s,X0,W,V,u)
%
%   INPUTS
%           A: system matrix A
%          Bu: system matrix Bu
%          Bw: system matrix Bw
%           r: process bias
%           C: system matrix C
%          Dv: system matrix Dv
%           s: measurement bias
%          X0: zonotope X0
%           W: zonotope W
%           V: zonotope V
%           u: input sequence (sequence along columns)
%              
%   OUTPUTS
%         Phi: reachable state sets expressed as a zonotope
%         Psi: reachable output sets expressed as a zonotope
%   Tildevars: matrices and vectors of the augmented system

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

nof_iterations = size(u,2);
nof_states = size(A,1);
nof_inputs = size(Bu,2);
nof_disturbances = size(Bw,2);


% Input sequence
utilde = reshape(u,nof_iterations*nof_inputs,1);


% Reachable state set variables
rtilde = cell(nof_iterations,1);
Atilde = cell(nof_iterations,1);
Butilde = cell(nof_iterations,nof_iterations);
Bwtilde = cell(nof_iterations,nof_iterations);
Buaux = cell(nof_iterations,1);
Bwaux = cell(nof_iterations,1);

rtilde{1} = r;
Atilde{1} = A;
Butilde{1,1} = Bu;
Bwtilde{1,1} = Bw;
for m=nof_iterations:-1:2
    Butilde{1,m} = zeros(nof_states,nof_inputs);
    Bwtilde{1,m} = zeros(nof_states,nof_disturbances);
end 
Buaux{1} = [Butilde{1,:}];
Bwaux{1} = [Bwtilde{1,:}];


for j=2:nof_iterations
       
    % Prediction model variables
    
    rtilde{j} = A*rtilde{j-1} + r;
    Atilde{j} = A*Atilde{j-1};
    Butilde{j,j} = Bu;
    Bwtilde{j,j} = Bw;
    for m=nof_iterations:-1:j+1
        Butilde{j,m} = zeros(nof_states,nof_inputs);
        Bwtilde{j,m} = zeros(nof_states,nof_disturbances);
    end    
    for m=j-1:-1:1
        Butilde{j,m} = A*Butilde{j-1,m};
        Bwtilde{j,m} = A*Bwtilde{j-1,m};
    end
    Buaux{j} = [Butilde{j,:}];
    Bwaux{j} = [Bwtilde{j,:}];
    
end

rtilde = vertcat(rtilde{:});
Atilde = vertcat(Atilde{:});
Butilde = vertcat(Buaux{:});
Bwtilde = vertcat(Bwaux{:});
Wtilde = repmat(W,nof_iterations);



% Reachable output set variables

stilde = repmat(s,nof_iterations,1);
Ctildeaux = repmat({C},nof_iterations,1);
Ctilde = blkdiag(Ctildeaux{:});
Dvtildeaux = repmat({Dv},nof_iterations,1);
Dvtilde = blkdiag(Dvtildeaux{:});
Vtilde = repmat(V,nof_iterations);


if(nargout==3)
    % Reachable sets
    Phi = Atilde*X0 + Bwtilde*Wtilde;
    Phi = Phi + Butilde*utilde;

    Psi = Ctilde*Phi + Dvtilde*Vtilde;
    Psi = Psi + stilde;
end

% Augmented variables
Tildevars.r = rtilde;
Tildevars.A = Atilde;
Tildevars.Bu = Butilde;
Tildevars.Bw = Bwtilde;
Tildevars.s = stilde;
Tildevars.C = Ctilde;
Tildevars.Dv = Dvtilde;
Tildevars.W = Wtilde;
Tildevars.V = Vtilde;

end



