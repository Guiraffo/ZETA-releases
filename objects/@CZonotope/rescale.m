function [Znew,csimid,csirad] = rescale(Z,method)
%RESCALE performs rescale over a constrained zonotope 
%
%   SYNTAX: Znew = RESCALE(Z)
%           Znew = RESCALE(Z,method)
%
%   INPUTS
%           Z: constrained zonotope as a structure object
%      method: 'LP': linear programming
%              'IA': interval arithmetic
%              (default is 'IA')
%   OUTPUT
%        Znew: the new constrained zonotope
%      csimid: the center of the new constrained box
%      csirad: the radius of the new constrained box

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

if(nargin<1)
    error('Not enough input arguments.');
elseif(nargin==1)
    method = 'IA';
end



if(strcmp(method, 'LP'))
   

    nof_generators = size(Z.G,2);
    nof_constraints = size(Z.A,1);

    
    csilower = -ones(nof_generators,1);
    csiupper = ones(nof_generators,1);
    
    OPTIONS.Display = 'off'; 
    
    for j=1:nof_generators

        f = zeros(nof_generators+1,1);
        f(j+1) = 1;
        A = [                      1,  zeros(1,nof_generators);
             -ones(nof_generators,1),      eye(nof_generators);
             -ones(nof_generators,1),     -eye(nof_generators)];
        b = [1; zeros(2*nof_generators,1)];
        Aeq = [zeros(nof_constraints,1), Z.A];
        beq = Z.b;

       
        x = optim.solvelp(f,A,b,Aeq,beq,[],[],OPTIONS);  csilower(j) =  x(j+1); % Revision 18-10-2024
        x = optim.solvelp(-f,A,b,Aeq,beq,[],[],OPTIONS); csiupper(j) =  x(j+1); % Revision 18-10-2024        
        
    end
    
    csimid = (csilower+csiupper)/2;
    csirad = (csiupper-csilower)/2;

    Znew_.c = Z.c + Z.G*csimid;
    Znew_.G = Z.G*diag(csirad);            
    Znew_.A = Z.A*diag(csirad);
    Znew_.b = Z.b - Z.A*csimid;

    Znew = CZonotope(Znew_.c,Znew_.G,Znew_.A,Znew_.b);
  
    
elseif(strcmp(method, 'IA'))
    
    [Znew_,~,~,~,csimid,csirad] = libCZon.Scale(libCZon.Wrapto(Z));
    Znew  = libCZon.Wrapfrom(Znew_);
    
    
else    
    
    error('Selected method for rescaling must be "LP" or "IA".');    

end


end





