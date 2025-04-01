function [Xhat,Xbar,cputimes] = czon_descriptor_estimator(dtsys,Xhat_prev,W,V,u,uprev,y,Xadm,OPTIONS)
%CZON_DESCRIPTOR_ESTIMATOR implements a state estimator for linear
%                          discrete-time systems using constrained
%                          zonotopes                       

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

% Setup algorithm flags based on input string
flags.pred = any(OPTIONS.alg=='p');
flags.updt = any(OPTIONS.alg=='u');

if(~(flags.pred || flags.updt))
    error('Invalid estimator operation mode.');
end

% Prediction step
tic;
if(flags.pred)
    Xbar = prediction(dtsys,Xhat_prev,W,u,uprev,Xadm);
else
    Xbar = Xhat_prev;
end
cputimes.pred = toc;

% Update step
tic;
if(flags.updt)
    Xhat = update(dtsys,Xbar,V,y);
else
    Xhat = Xbar;
end
cputimes.updt = toc;

% Complexity reduction
tic;
Xhat = reduction(Xhat,OPTIONS.ngmax,OPTIONS.ncmax);
cputimes.redu = toc;

% Total computational time
cputimes.total = cputimes.pred + cputimes.updt + cputimes.redu;

end

function Xbar = prediction(dtsys,X,W,u,uprev,Xadm)
% Prediction step of descriptor systems using CZs

SVDvars = dtsys.SVDvars;
Tinv = inv(SVDvars.T);

% Move the admissible set to the transformed space
Xadm_SVD = Tinv*Xadm;
cadm_SVD = Xadm_SVD.c(SVDvars.nof_nonzeros+1:end,1);
Gadm_SVD = Xadm_SVD.G(SVDvars.nof_nonzeros+1:end,:);

% Move to the transformed space
X_SVD = Tinv*X;

% Prediction in the transformed space
Xbar_SVD =   CZonotope([SVDvars.Atilde*X_SVD.c + SVDvars.Btilde*uprev + SVDvars.Bwtilde*W.c;
                                                                                   cadm_SVD],...
                         [SVDvars.Atilde*X_SVD.G, SVDvars.Bwtilde*W.G, zeros(size(SVDvars.Atilde,1),size(Gadm_SVD,2)), zeros(size(SVDvars.Atilde,1),size(W.G,2));
                          zeros(size(Gadm_SVD,1),size(X_SVD.G,2)), zeros(size(Gadm_SVD,1),size(W.G,2)), Gadm_SVD, zeros(size(Gadm_SVD,1),size(W.G,2))],...
                         [                          X_SVD.A, zeros(size(X_SVD.A,1),size(W.A,2)), zeros(size(X_SVD.A,1),size(Gadm_SVD,2)), zeros(size(X_SVD.A,1),size(W.G,2));
                          zeros(size(W.A,1),size(X_SVD.A,2)),                               W.A,     zeros(size(W.A,1),size(Gadm_SVD,2)), zeros(size(W.A,1),size(W.G,2));
                          SVDvars.Ahat*[SVDvars.Atilde*X_SVD.G, SVDvars.Bwtilde*W.G, zeros(size(SVDvars.Atilde,1),size(Gadm_SVD,2));
                          zeros(size(Gadm_SVD,1),size(X_SVD.G,2)), zeros(size(Gadm_SVD,1),size(W.G,2)), Gadm_SVD], SVDvars.Bwhat*W.G],...
                          [X_SVD.b;
                               W.b;
               -SVDvars.Ahat*[SVDvars.Atilde*X_SVD.c + SVDvars.Btilde*uprev + SVDvars.Bwtilde*W.c;
                                                   cadm_SVD] - SVDvars.Bhat*u - SVDvars.Bwhat*W.c]);
                                               
Xbar = SVDvars.T*Xbar_SVD; % Return to the original space                                                  

end

function Xhat = update(dtsys,Xbar,V,y)
% Update step of descriptor systems using CZs

Y = y - dtsys.Dv*V;
Xhat = intersection(Xbar,Y,dtsys.C);  

end



