function [e_minus,e_plus] = get_error_bounds(f,vert,hz,OPTIONS,lambda)
%GET_ERROR_BOUNDS generates the linearization error bounds using DC
%                 programming principles. Internal use only.

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

nof_vertices = size(vert,2);

err_minor = [];
err_major = [];

for j=1:nof_vertices

    z = vert(:,j);
    err_minor(:,j) = error_minorant(f,z,hz,OPTIONS,lambda);
    err_major(:,j) = error_majorant(f,z,hz,OPTIONS,lambda);

end

e_minus = min(err_minor,[],2);
e_plus = max(err_major,[],2);

end

%% DC programming: Error minorant and majorant
function out = error_minorant(f,z,hz,OPTIONS,lambda)

    linearized_f_at_xw = linearized_function(f,z,hz);  
    switch OPTIONS.decomposition
        case 'exact'        
            fconvxb_at_xw = f.dcB(z);
            linearized_fconvxa_at_xw = linearized_function_convxa(f,z,hz);
        case 'aBB'            
            fconvxb_at_xw = dcprog.evalfunction_convxbABB(f,z,lambda);
            linearized_fconvxa_at_xw = linearized_function_convxaABB(f,z,hz,lambda);
        otherwise
            error('Invalid decomposition mode.');
    end

    out = linearized_fconvxa_at_xw - fconvxb_at_xw - linearized_f_at_xw;
    
end
function out = error_majorant(f,z,hz,OPTIONS,lambda)

    linearized_f_at_xw = linearized_function(f,z,hz);
    switch OPTIONS.decomposition
        case 'exact'     
            fconvxa_at_xw = f.dcA(z);
            linearized_fconvxb_at_xw = linearized_function_convxb(f,z,hz);
        case 'aBB'
            fconvxa_at_xw = dcprog.evalfunction_convxaABB(f,z,lambda);
            linearized_fconvxb_at_xw = linearized_function_convxbABB(f,z,hz,lambda);  
        otherwise
            error('Invalid decomposition mode.');
    end            
    
    out = fconvxa_at_xw - linearized_fconvxb_at_xw - linearized_f_at_xw;
    
end
%% DC programming: Evaluation of linearized functions f, fa, fb
function out = linearized_function(f,z,hz)
    out = f.eval(hz) + f.Jacob(hz)*(z-hz);
end
function out = linearized_function_convxa(f,z,hz)
    out = f.dcA(hz) + f.dcA_Jacob(hz)*(z-hz);
end
function out = linearized_function_convxb(f,z,hz)
    out = f.dcB(hz) + f.dcB_Jacob(hz)*(z-hz);
end
%% DC programming: Evaluation of fa and fb, and respective Jacobians, using the alpha-BB method, plus linearized versions
% function out = evalfunction_convxaABB(f,z,lambda)
% % fconvxa for system DC decomposition using alphaBB method: fa = f + fb, fb = 0.5*lambda*(z.'z)
% 
% out = f.eval(z) + 0.5*lambda*(z.'*z);
% 
% end
% function out = evalfunction_convxaABB_Jacob(f,z,lambda)
% % fconvxa (Jacobian) for system DC decomposition using alphaBB method: fa = f + fb, fb = 0.5*lambda*(z.'z)
% 
% out = f.Jacob(z) + lambda*z.';
% 
% end
% function out = evalfunction_convxbABB(f,z,lambda)
% % fconvxb for system DC decomposition using alphaBB convexification: fa = f + fb, fb = 0.5*lambda*(z.'z)
% 
% out = 0.5*lambda*(z.'*z);
% 
% end
% function out = evalfunction_convxbABB_Jacob(f,z,lambda)
% % fconvxb (Jacobian) for system DC decomposition using alphaBB method: fa = f + fb, fb = 0.5*lambda*(z.'z)
% 
% out = lambda*z.';
% 
% end
function out = linearized_function_convxaABB(f,z,hz,lambda)
    out = dcprog.evalfunction_convxaABB(f,hz,lambda) + dcprog.evalfunction_convxaABB_Jacob(f,hz,lambda)*(z-hz);
end
function out = linearized_function_convxbABB(f,z,hz,lambda)
    out = dcprog.evalfunction_convxbABB(f,hz,lambda) + dcprog.evalfunction_convxbABB_Jacob(f,hz,lambda)*(z-hz);
end
