function [Xhat,Xbar,cputimes] = zon_meanvalue_estimator(dtsys,Xhat_prev,W,V,u,y,OPTIONS)
%ZON_MEANVALUE_ESTIMATOR implements a state estimator for nonlinear
%                        discrete-time systems using zonotopes based on the
%                        mean value extension

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
    Xbar = prediction(dtsys.f,Xhat_prev,W,u,dtsys.Ts);
else
    Xbar = Xhat_prev;
end
cputimes.pred = toc;

% Update step
tic
if(flags.updt)
    Xhat = update(dtsys.g,Xbar,V,y);
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

function Xbar = prediction(f,X,W,u,Ts)

% Implementation of Theorem 4 from "Guaranteed State Estimation by Zonotopes" (Alamo et al, 2005)

% ----------------------------- Interval matrices M --------------------------------

Xhull = intervalhull(X);

if(isempty(W))
    W = Zonotope;
    Whull = [];
else
    Whull = intervalhull(W);
end

%Jacobians
Jacobf_x = f.Jacob(Xhull,Whull,u,Ts);
Jacobf_w = f.Jacob(X.c,Whull,u,Ts);
Jacobf_x = Jacobf_x(:,1:size(Xhull,1));
Jacobf_w = Jacobf_w(:,size(Xhull,1)+1:end);
           
M_x = interval(Jacobf_x*X.G); % Interval call just in case the result is not an interval, so mid and diam still work
M_w = interval(Jacobf_w*W.G); % Interval call just in case the result is not an interval, so mid and diam still work

% -----------------------------------------------------------------------

% --------------------------- Zonotope Xbar ------------------------------


Zinc_x = Zonotope.inclusion(zeros(size(M_x,1),1), M_x);
Zinc_w = Zonotope.inclusion(zeros(size(M_w,1),1), M_w);

Xbar = f.eval(X.c,W.c,u,Ts) + Zinc_x + Zinc_w;


% -----------------------------------------------------------------------

end

function Xhat = update(g,Xbar,V,y)
% Nonlinear measurements by Alamo et al (2005) and intersection by Bravo et al (2006)

Xtilde = Xbar;
if(isempty(V))
    V = Zonotope;
    Vhull = [];
else
    Vhull = intervalhull(V);
end

% For loop

for i=1:length(y)
    
    % Compute the strip of consistent states
    
    % Get interval Jacobian and the strip parameter p
    Xhull = intervalhull(Xtilde);
    Jacobg = g.Jacob(Xhull,Vhull);
    Jacobg_x = Jacobg(:,1:size(Xhull,1));

    p = mid(Jacobg_x(i,:).');

    % Computation of the strip parameters d and sigma

    if(sum(diam(Jacobg(i,:)))==0) % If the Jacobian is not an interval, it means this measurement is linear

        Dvi = mid(Jacobg(i,size(Xhull,1)+1:end));
        DvVi = Dvi.'*V;
        d = y(i) - DvVi.c(i);
        sigma = sum(rad(intervalhull(DvVi)));

    else % Otherwise, use the general formula in Alamo et al (2005)

        outputhull = g.eval(Xhull,Vhull);
        stripinterval = p.'*Xhull - outputhull(i);
        d = y(i) + mid(stripinterval);
        sigma = rad(stripinterval);

    end
    
    % Intersect with the strip of consistent states
    Xtilde = intersection(Xtilde,strip(p, d, sigma));
    
end
    
Xhat = Xtilde;

end

