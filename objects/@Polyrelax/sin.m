function Z = sin(X)
%SIN returns the sine of a polyrelax object
%
%   SYNTAX: Z = SIN(X)
%
%   INPUTS
%           X: the input polyrelax object
%
%   OUTPUT
%           Z: the resulting polyrelax object

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

X = Polyrelax(X);
Z = Polyrelax(sin(X.x));

xinfx1 = 2*pi;
xinfx2 = 3*pi;
xmin1 = 3*pi/2;
xmin2 = 7*pi/2;

xL = inf(X.x);
xU = sup(X.x);

% % Check if input is not contained in [3*pi/2,7*pi/2]
% if((xL<xmin1)||(xL>xmin2)||(xU<xmin1)||(xU>xmin2))
%     error(' The sine function is currently implemented only for the interval [3*pi/2,7*pi/2].');
% end


% Get object indexes
Xind = X.i;
Zind = Z.i;

% Get convex polyhedral relaxation
[Hconvx,kconvx] = convexrelax(xL,xU,xmin1,xmin2,xinfx1,xinfx2,Xind,Zind);

% Get concave polyhedral relaxation
[Hconcv,kconcv] = convexrelax(-xU,-xL,xmin1,xmin2,xinfx1,xinfx2,Xind,Zind);
Hconcv(:,Zind) = -Hconcv(:,Zind);
Hconcv(:,Xind) = -Hconcv(:,Xind);

% Total polyhedral enclosure    
Hadd = [Hconvx;Hconcv];
kadd = [kconvx;kconcv];

% Adds the inequality constraints to polyrelax Hrep
Polyrelax.addHrepHrow(Hadd,kadd);

% %% Debugging section
%        
% keyboard;
% xsample = xL:0.01:xU;
% zsample = sin(xsample);
% zL = inf(Z.x);
% zU = sup(Z.x);
% %P = Polyhedron('A',Hconvx,'b',kconvx,'LB',[xL;-Inf],'UB',[xU;zU]);
% P = Polyhedron('A',Hadd,'b',kadd);
% figure
% plot(xsample,zsample,'LineWidth',2);
% hold on
% plot(P,'Color','b','Alpha',0.1);
   

end

% Convex relaxation function
function [Hconvx,kconvx] = convexrelax(xL,xU,xmin1,xmin2,xinfx1,xinfx2,Xind,Zind)

% Compute auxiliary variables
[n1,n2] = n12(xL,xU);
[wL,wU,vL,vU] = getWV(xL,xU,n1,n2,xmin1,xmin2);


% Verify the existence of the three intervals in pag 120 of Scott's thesis
flags.interval1 = (xL - 2*(n1 - 1)*pi <= xmin2);
flags.interval2 = (xU - 2*(n2 - 1)*pi >= xmin1);
flags.otherwise = ((xU - 2*(n1 - 1)*pi > xmin2)||(xL - 2*(n2 - 1)*pi < xmin1));


% If 'otherwise' interval exists
if(flags.otherwise)  
    % z > -1
    Hother = zeros(1,Zind);
    kother = zeros(1,1);               
    Hother(1,Zind) = -1;
    kother(1,1) = 1;    
else % If it does not exist    
    Hother = zeros(0,Zind);
    kother = zeros(0,1);   
end

% If interval 1 exists
if(flags.interval1)
    [Hinter1,kinter1] = etarelax(wL,wU,xmin1,xmin2,xinfx1,xinfx2,Xind,Zind,-2*(n1-1)*pi);
else
    Hinter1 = zeros(0,Zind);
    kinter1 = zeros(0,1);    
end

% If interval 2 exists
if(flags.interval1)
    [Hinter2,kinter2] = etarelax(vL,vU,xmin1,xmin2,xinfx1,xinfx2,Xind,Zind,-2*(n2-1)*pi);
else
    Hinter2 = zeros(0,Zind);
    kinter2 = zeros(0,1);    
end

Hconvx = [Hother;Hinter1;Hinter2];
kconvx = [kother;kinter1;kinter2];


end

% Eta relaxation function
function [H,k] = etarelax(xL,xU,xmin1,xmin2,xinfx1,xinfx2,Xind,Zind,xbias)
% Implements a polyhedral version of the eta convex relaxation function from p. 119 in Scott's thesis
% Inputs: [xL,xU]: evaluated interval
%         xmin1,xmin2: argmin of the sin function in [3pi/2, 7pi/2]
%         xinfx1, xinfx2: inflection points of the sin function in [3pi/2, 7pi/2]
%         Xind = index of polyrelax object X (input object)
%         Zind = index of polyrelax object Z (output object)
%         xbias = constant that translates the argument x (x + xbias)
% Output: [H,k]: associated Qcv

% Compute auxiliary variables

% Get x prime variables
persistent OPTIONS;
if(isempty(OPTIONS))
    %OPTIONS = optimoptions('fsolve');
    %OPTIONS = optimoptions('fsolve','Display','off');
    %OPTIONS = optimset('Display','off');
    OPTIONS = optimset('fzero');
    OPTIONS.MyTOL = 1e-12;
end
% x_pL = fsolve(@(x) trigconstrU(x,xU), 0.5*(xmin1+xinfx1), OPTIONS);
% x_pU = fsolve(@(x) trigconstrL(x,xL), 0.5*(xmin2+xinfx2), OPTIONS);

% Get x star variables
if(xU<=xinfx1)
    x_sL = xinfx1;
else
    %x_sL = x_pL;
    %x_sL = fsolve(@(x) trigconstrU(x,xU), 0.5*(xmin1+xinfx1), OPTIONS);
    x_sL = fzero(@(x) trigconstrU(x,xU), [xmin1-OPTIONS.MyTOL,xinfx1+OPTIONS.MyTOL], OPTIONS);
end
if(xL>=xinfx2)
    x_sU = xinfx2;
else
    %x_sU = x_pU;
    %x_sU = fsolve(@(x) trigconstrL(x,xL), 0.5*(xmin2+xinfx2), OPTIONS);
    x_sU = fzero(@(x) trigconstrL(x,xL), [xinfx2-OPTIONS.MyTOL,xmin2+OPTIONS.MyTOL], OPTIONS);
end

% Get x1 and x2 variables
x1 = mid3([xL,xU,x_sL]);
x2 = mid3([xL,xU,x_sU]);

%{
There are 3 intervals to be analysed: left, middle, right. 
The function is convex in the left interval.
The function is convex in the right interval.
The secant must be obtained for the middle interval.
In some cases, not all the 3 exist.
Need to verify if they exist for the input interval.
%}

% Flags for the existence of the intervals
flags.left = 1;
flags.right = 1;
flags.middle = 1;
if(x1==xU) % Only the left interval exists
    flags.middle = 0;
    flags.right = 0;
elseif(x2==xL) % Only the right interval exists
    flags.left = 0;
    flags.middle = 0;
else
    if(x1==xL) % Only the middle and right intervals exist
        flags.left = 0;
    end
    if(x2==xU) % Only the left and middle intervals exist
        flags.right = 0;
    end
end

% Obtain the secant-based Hrep for the middle interval
if(flags.middle) % If the middle interval exists
    Hsec = zeros(1,Zind);
    ksec = zeros(1,1);

    % Convex secant
    % z > (sin(x2) - sin(x1))/(x2 - x1))*(x - x1) + sin(x1)
    secant_slope = (sin(x2)-sin(x1))/(x2-x1);
    Hsec(1,Xind) = secant_slope;
    Hsec(1,Zind) = -1; 
    ksec(1,1) = secant_slope*x1 - sin(x1);  

else % If it does not exist
    Hsec = zeros(0,Zind);
    ksec = zeros(0,1);
end

% Linear approximation for left interval
if(flags.left) % If the left interval exists (uses linearization)
    HLef = zeros(3,Zind);
    kLef = zeros(3,1);               
    xLefM = 0.5*(xL+x1);

    % Linearized convex lowerbound for Right interval
    % z > diff(sin(x))|(x=xL)*(x - xL) + sin(xL)    
    % z > diff(sin(x))|(x=xmid)*(x - xmid) + sin(xmid)
    % z > diff(sin(x))|(x=x1)*(x - x1) + sin(x1)
    HLef(1,Xind) = cos(xL);    
    HLef(2,Xind) = cos(xLefM);
    HLef(3,Xind) = cos(x1);    
    HLef(1,Zind) = -1;
    HLef(2,Zind) = -1;
    HLef(3,Zind) = -1;    
    kLef(1,1) = cos(xL)*xL - sin(xL);
    kLef(2,1) = cos(xLefM)*xLefM - sin(xLefM);
    kLef(3,1) = cos(x1)*x1 - sin(x1);    
   
else % If it does not exist    
    HLef = zeros(0,Zind);
    kLef = zeros(0,1);   
end

% Linear approximation for right interval
if(flags.right) % If the right interval exists (uses linearization)
    HRig = zeros(3,Zind);
    kRig = zeros(3,1);               
    xRigM = 0.5*(x2+xU);

    % Linearized convex lowerbound for Right interval
    % z > diff(sin(x))|(x=xmid)*(x - xmid) + sin(xmid)    
    % z > diff(sin(x))|(x=xU)*(x - xU) + sin(xU)
    % z > diff(sin(x))|(x=x2)*(x - x2) + sin(x2)    
    HRig(1,Xind) = cos(xRigM);    
    HRig(2,Xind) = cos(xU);
    HRig(3,Xind) = cos(x2);    
    HRig(1,Zind) = -1;
    HRig(2,Zind) = -1;
    HRig(3,Zind) = -1;    
    kRig(1,1) = cos(xRigM)*xRigM - sin(xRigM);
    kRig(2,1) = cos(xU)*xU - sin(xU);
    kRig(3,1) = cos(x2)*x2 - sin(x2);    
   
else % If it does not exist    
    HRig = zeros(0,Zind);
    kRig = zeros(0,1);   
end
 
H = [Hsec;HLef;HRig];
k = [ksec;kLef;kRig];

% Incorporates xbias
bias = zeros(Zind,1);
bias(Xind,1) = xbias;
k = k - H*bias;

end

% Auxiliary functions
function [wL,wU,vL,vU] = getWV(xL,xU,n1,n2,xmin1,xmin2)
% Computes the intervals Z and Y (in this code, W and V, respectively)
wL = xL - 2*(n1-1)*pi;
wU = min([xU - 2*(n1-1)*pi, xmin2]);
vL = xmin1;
vU = xU - 2*(n2-1)*pi;
end
function [n1,n2] = n12(xL,xU)
% Obtain n1 and n2 using the n(x) function
n1 = floor(nx(xL));
nxU = nx(xU);
if(mod(nxU,1)==0) % If its an integer
    n2 = nxU-1;
else
    n2 = floor(nxU);
end
end
function out = nx(x)
% Evaluates the n(x) function
out = x/(2*pi) + (1/4);
end
function out = mid3(x)
% Returns the mid point of a 3 elements array
x = sort(x);
out = x(2);
end
function out = trigconstrU(x,xU)
% Implements (xU - x)*cos(x) + sin(x) - sin(xU) = 0
out = (xU - x)*cos(x) + sin(x) - sin(xU);
end
function out = trigconstrL(x,xL)
% Implements (x - xL)*cos(x) - sin(x) + sin(xL) = 0
out = (x - xL)*cos(x) - sin(x) + sin(xL);
end