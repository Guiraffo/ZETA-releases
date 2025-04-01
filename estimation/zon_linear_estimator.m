function [Xhat,Xbar,cputimes] = zon_linear_estimator(dtsys,Xhat_prev,W,V,u,y,OPTIONS)
%ZON_LINEAR_ESTIMATOR implements a state estimator for linear discrete-time
%                     systems using zonotopes                       

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
    Xbar = prediction(dtsys,Xhat_prev,W,u);
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
Xhat = reduction(Xhat,OPTIONS.ngmax);
cputimes.redu = toc;

% Total computational time
cputimes.total = cputimes.pred + cputimes.updt + cputimes.redu;

end

function Xbar = prediction(dtsys,Xhat_prev,W,u)
% Linear prediction step of linear systems using zonotopes

Xbar = dtsys.A*Xhat_prev + dtsys.Bu*u + dtsys.Bw*W;

end

function Xhat = update(dtsys,Xbar,V,y)
% Linear update step of linear systems using zonotopes based on Bravo et
% al. (2006)

Xtilde = Xbar;
DvV = dtsys.Dv*V;
DvVhull = intervalhull(DvV);  

% Iterative intersection with strips

for i=1:length(y)
    
    % Compute the parameters of the strip of consistent states
    p = dtsys.C(i,:).';
    d = y(i) - DvV.c(i);
    sigma = rad(DvVhull(i,:));
    
    Xtilde = intersection(Xtilde,Strip(p, d, sigma));
    
end
       
Xhat = Xtilde;

end



