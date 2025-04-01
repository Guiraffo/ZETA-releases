function Zimpr = partopeimprove(Z,f,P_vert,hz,OPTIONS,lambda)
%PARTOPEIMPROVE implements the parallelotope improvement step in the
%               zonotope-based DC programming enclosure propagation
%               (Alamo et al., 2008). Internal use only.

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

Phat = partopebound(Z);

E = inv(Phat.G); % Easier way to get the half-spaces matrix

[p_minus,p_plus] = get_improvement_bounds(f,E,P_vert,hz,OPTIONS,lambda);
ImprBound = Zonotope(0.5*(p_plus+p_minus), 0.5*diag(p_plus-p_minus));


% Iterative intersection with strips
Zimpr = Z;
for i=1:size(E,1)
    
    % Compute the parameters of the associated strip 
    p = E(i,:).';
    d = ImprBound.c(i);
    sigma = sum(abs(ImprBound.G(i,:)));
    
    Zimpr = intersection(Zimpr,Strip(p, d, sigma));
    
end

end

function [p_minus,p_plus] = get_improvement_bounds(f,E,vert,hz,OPTIONS,lambda)

    nof_vertices = size(vert,2);

    improvement_minor = [];
    improvement_major = [];

    for j=1:nof_vertices
        
        z = vert(:,j);
        improvement_minor(:,j) = improvement_minorant(f,E,z,hz,OPTIONS,lambda);
        improvement_major(:,j) = improvement_majorant(f,E,z,hz,OPTIONS,lambda);

    end

    p_minus = min(improvement_minor,[],2);
    p_plus  = max(improvement_major,[],2);

end

function out = improvement_minorant(f,E,z,hz,OPTIONS,lambda)

    fconvxbE_at_xw = evalfunction_convxbE(f,E,z,OPTIONS,lambda);
    linearized_fconvxaE_at_xw = linearized_function_convxaE(f,E,z,hz,OPTIONS,lambda);

    out = linearized_fconvxaE_at_xw - fconvxbE_at_xw;
    
end
function out = improvement_majorant(f,E,z,hz,OPTIONS,lambda)

    fconvxaE_at_xw = evalfunction_convxaE(f,E,z,OPTIONS,lambda);
    linearized_fconvxbE_at_xw = linearized_function_convxbE(f,E,z,hz,OPTIONS,lambda);
    
    out = fconvxaE_at_xw - linearized_fconvxbE_at_xw;
    
end
function out = evalfunction_convxaE(f,E,z,OPTIONS,lambda)

    switch OPTIONS.decomposition
        case 'exact' 
            convxa = f.dcA(z);
            convxb = f.dcB(z);
        case 'aBB'
            convxa = dcprog.evalfunction_convxaABB(f,z,lambda);
            convxb = dcprog.evalfunction_convxbABB(f,z,lambda);
        otherwise
            error('Invalid decomposition mode.');
    end 
    
    f_size = size(E,2);    
    out = zeros(f_size,1);
    for i=1:f_size
        for j=1:f_size
            if(E(i,j)>=0)
                out(i,:) = out(i,:) + E(i,j)*convxa(j,1);
            else
                out(i,:) = out(i,:) - E(i,j)*convxb(j,1);
            end
        end
    end
end
function out = evalfunction_convxbE(f,E,z,OPTIONS,lambda)

    switch OPTIONS.decomposition
        case 'exact' 
            convxa = f.dcA(z);
            convxb = f.dcB(z);
        case 'aBB'
            convxa = dcprog.evalfunction_convxaABB(f,z,lambda);
            convxb = dcprog.evalfunction_convxbABB(f,z,lambda);
        otherwise
            error('Invalid decomposition mode.');
    end 
    
    f_size = size(E,2);    
    out = zeros(f_size,1);
    for i=1:f_size
        for j=1:f_size
            if(E(i,j)>=0)
                out(i,:) = out(i,:) + E(i,j)*convxb(j,1);
            else
                out(i,:) = out(i,:) - E(i,j)*convxa(j,1);
            end
        end
    end
end
function out = linearized_function_convxaE(f,E,z,hz,OPTIONS,lambda)

    switch OPTIONS.decomposition
        case 'exact' 
            fconvxa_Jacob = f.dcA_Jacob(hz);
            fconvxb_Jacob = f.dcB_Jacob(hz);
        case 'aBB'
            fconvxa_Jacob = dcprog.evalfunction_convxaABB_Jacob(f,hz,lambda);
            fconvxb_Jacob = dcprog.evalfunction_convxbABB_Jacob(f,hz,lambda);
        otherwise
            error('Invalid decomposition mode.');
    end 
    
    fconvxaE_Jacob = zeros(size(E,1),size(fconvxa_Jacob,2));
    for i=1:size(E,1)
        for k=1:size(fconvxa_Jacob,2)
            for j=1:size(fconvxa_Jacob,1)
                if(E(i,j)>=0)
                    fconvxaE_Jacob(i,k) = fconvxaE_Jacob(i,k) + E(i,j)*fconvxa_Jacob(j,k);
                else
                    fconvxaE_Jacob(i,k) = fconvxaE_Jacob(i,k) - E(i,j)*fconvxb_Jacob(j,k);
                end            
            end
        end
    end
    
    out = evalfunction_convxaE(f,E,hz,OPTIONS,lambda) + fconvxaE_Jacob*(z-hz);
end
function out = linearized_function_convxbE(f,E,z,hz,OPTIONS,lambda)

    switch OPTIONS.decomposition
        case 'exact' 
            fconvxa_Jacob = f.dcA_Jacob(hz);
            fconvxb_Jacob = f.dcB_Jacob(hz);
        case 'aBB'
            fconvxa_Jacob = dcprog.evalfunction_convxaABB_Jacob(f,hz,lambda);
            fconvxb_Jacob = dcprog.evalfunction_convxbABB_Jacob(f,hz,lambda);
        otherwise
            error('Invalid decomposition mode.');
    end 
    
    fconvxbE_Jacob = zeros(size(E,1),size(fconvxa_Jacob,2));
    for i=1:size(E,1)
        for k=1:size(fconvxa_Jacob,2)
            for j=1:size(fconvxa_Jacob,1)
                if(E(i,j)>=0)
                    fconvxbE_Jacob(i,k) = fconvxbE_Jacob(i,k) + E(i,j)*fconvxb_Jacob(j,k);
                else
                    fconvxbE_Jacob(i,k) = fconvxbE_Jacob(i,k) - E(i,j)*fconvxa_Jacob(j,k);
                end            
            end
        end
    end
    
    out = evalfunction_convxbE(f,E,hz,OPTIONS,lambda) + fconvxbE_Jacob*(z-hz);
end

