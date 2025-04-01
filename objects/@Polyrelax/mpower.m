function Z = mpower(X,a)
%^ or MPOWER returns the power of a Polyrelax object
%
%   SYNTAX: Z = X^a
%           Z = MPOWER(X,a)
%
%   INPUTS
%           Z: the Polyrelax object
%           a: the exponent (a positive number)
%
%   OUTPUT
%           Z: the resulting Polyrelax object

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

% Input check
if(~isnumeric(a))
    error('The second argument must be a numeric variable.');
end
if(numel(a)~=1)
    error('The second argument must be a scalar.');
end
if(a < 0)
    error('Only positive exponents are currently implemented.');
end


X = Polyrelax(X);
Z = Polyrelax(X.x^a);


% Adds rows to Polyrelax Hrep
if(rem(a, 2)==0) % Even exponent
    
power_even(X,Z,a);
    
elseif(rem(a, 2)==1) % Odd exponent

power_odd(X,Z,a);    
    
else
    disp('WARNING: positive non-integer exponents are experimental.');
power_noninteger(X,Z,a);  
end


% 21-05-2024: added approximation modes for even integer power
    
end

function power_even(X,Z,a)

xL = inf(X.x);
xU = sup(X.x);

% Concave upperbound
% z < ((xU^a - xL^a)/(xU - xL))*(x - xL) + xL^a
Hconcv = zeros(1,Z.i);
kconcv = zeros(1,1);
concv_slope = ((xU^a - xL^a)/(xU - xL));
Hconcv(1,X.i) = -concv_slope;
Hconcv(1,Z.i) = 1; 
kconcv(1,1) = -concv_slope*xL + xL^a;

% Convex lowerbound
% z > diff(x^a,x)|(x=xL)*(x - xL) + xL^a
% z > diff(x^a,x)|(x=xmid)*(x - xmid) + xmid^a
% z > diff(x^a,x)|(x=xU)*(x - xU) + xU^a
switch toolsettings.polyrelax_approxmode
    case 0    
        Hconvx = zeros(3,Z.i);
        kconvx = zeros(3,1);
        xmid = 0.5*(xL+xU);    
        Hconvx(1,X.i) = a*xL^(a-1);
        Hconvx(2,X.i) = a*xmid^(a-1);
        Hconvx(3,X.i) = a*xU^(a-1);
        Hconvx(1,Z.i) = -1;
        Hconvx(2,Z.i) = -1;
        Hconvx(3,Z.i) = -1;
        kconvx(1,1) = a*xL^(a-1)*xL - xL^a;
        kconvx(2,1) = a*xmid^(a-1)*xmid - xmid^a;
        kconvx(3,1) = a*xU^(a-1)*xU - xU^a;
    case 2
        Hconvx = zeros(2,Z.i);
        kconvx = zeros(2,1);  
        xlow = (2/3)*xL + (1/3)*xU;
        xupp = (1/3)*xL + (2/3)*xU;
        Hconvx(1,X.i) = a*xlow^(a-1);
        Hconvx(2,X.i) = a*xupp^(a-1);
        Hconvx(1,Z.i) = -1;
        Hconvx(2,Z.i) = -1;
        kconvx(1,1) = a*xlow^(a-1)*xlow - xlow^a;
        kconvx(2,1) = a*xupp^(a-1)*xupp - xupp^a;        
    otherwise
        error('Invalid Polyrelax approximation mode.')
end        
        
    
% Adds the inequality constraints to Polyrelax Hrep
Polyrelax.addHrepHrow([Hconcv;Hconvx],[kconcv;kconvx]);

end

function power_odd(X,Z,a)

    xL = inf(X.x);
    xU = sup(X.x);   
    
    %linearapprox = 1;
    linearapprox = 0;
    
    
    
    if(xL >= 0) % Positive interval: convex function
        
        xmid = 0.5*(xL+xU);
        Hadd = zeros(4,Z.i);
        kadd = zeros(4,1);

        % Concave upperbound
        % z < ((xU^a - xL^a)/(xU - xL))*(x - xL) + xL^a
        concv_slope = ((xU^a - xL^a)/(xU - xL));
        Hadd(1,X.i) = -concv_slope;
        Hadd(1,Z.i) = 1; 
        kadd(1,1) = -concv_slope*xL + xL^a;

        % Convex lowerbound
        % z > diff(x^a,x)|(x=xL)*(x - xL) + xL^a
        % z > diff(x^a,x)|(x=xmid)*(x - xmid) + xmid^a
        % z > diff(x^a,x)|(x=xU)*(x - xU) + xU^a
        Hadd(2,X.i) = a*xL^(a-1);
        Hadd(3,X.i) = a*xmid^(a-1);
        Hadd(4,X.i) = a*xU^(a-1);
        Hadd(2,Z.i) = -1;
        Hadd(3,Z.i) = -1;
        Hadd(4,Z.i) = -1;
        kadd(2,1) = a*xL^(a-1)*xL - xL^a;
        kadd(3,1) = a*xmid^(a-1)*xmid - xmid^a;
        kadd(4,1) = a*xU^(a-1)*xU - xU^a;     
        
    elseif(xU <= 0) % Negative interval: concave function
        
        xmid = 0.5*(xL+xU);
        Hadd = zeros(4,Z.i);
        kadd = zeros(4,1);

        % Concave upperbound
        % z < diff(x^a,x)|(x=xL)*(x - xL) + xL^a
        % z < diff(x^a,x)|(x=xmid)*(x - xmid) + xmid^a
        % z < diff(x^a,x)|(x=xU)*(x - xU) + xU^a
        Hadd(1,X.i) = -a*xL^(a-1);
        Hadd(2,X.i) = -a*xmid^(a-1);
        Hadd(3,X.i) = -a*xU^(a-1);
        Hadd(1,Z.i) = 1;
        Hadd(2,Z.i) = 1;
        Hadd(3,Z.i) = 1;
        kadd(1,1) = -a*xL^(a-1)*xL + xL^a;
        kadd(2,1) = -a*xmid^(a-1)*xmid + xmid^a;
        kadd(3,1) = -a*xU^(a-1)*xU + xU^a;          
        
        % Convex lowerbound        
        % z > ((xU^a - xL^a)/(xU - xL))*(x - xL) + xL^a
        convx_slope = ((xU^a - xL^a)/(xU - xL));
        Hadd(4,X.i) = convx_slope;
        Hadd(4,Z.i) = -1; 
        kadd(4,1) = convx_slope*xL - xL^a;        
        
        
    else % Interval that contains zero
        
        if(linearapprox)
        
        % This is an experimental linear approximation that does not
        % require to find the roots of the polynomials
        
        % Splits the interval into 'left' and 'right intervals
        xLefL = xL;
        xLefU = 0;
        xRigL = 0;
        xRigU = xU;
        
        xLefmid = 0.5*(xLefL+xLefU);
        xRigmid = 0.5*(xRigL+xRigU);
    
        HLef = zeros(4,Z.i); 
        kLef = zeros(4,1);
        HRig = zeros(4,Z.i);
        kRig = zeros(4,1);

        
        % Envelope for the left interval
        % Concave upperbound
        % z < diff(x^a,x)|(x=xL)*(x - xL) + xL^a
        % z < diff(x^a,x)|(x=xmid)*(x - xmid) + xmid^a
        % z < diff(x^a,x)|(x=xU)*(x - xU) + xU^a
        HLef(1,X.i) = -a*xLefL^(a-1);
        HLef(2,X.i) = -a*xLefmid^(a-1);
        HLef(3,X.i) = -a*xLefU^(a-1);
        HLef(1,Z.i) = 1;
        HLef(2,Z.i) = 1;
        HLef(3,Z.i) = 1;
        kLef(1,1) = -a*xLefL^(a-1)*xLefL + xLefL^a;
        kLef(2,1) = -a*xLefmid^(a-1)*xLefmid + xLefmid^a;
        kLef(3,1) = -a*xLefU^(a-1)*xLefU + xLefU^a;          
        
        % Convex lowerbound        
        % z > ((xU^a - xL^a)/(xU - xL))*(x - xL) + xL^a
        concv_slope_Lef = ((xLefU^a - xLefL^a)/(xLefU - xLefL));
        HLef(4,X.i) = concv_slope_Lef;
        HLef(4,Z.i) = -1; 
        kLef(4,1) = concv_slope_Lef*xLefL - xLefL^a;           
        
            
        % Envelope for the right interval
        % Concave upperbound
        % z < ((xU^a - xL^a)/(xU - xL))*(x - xL) + xL^a
        concv_slope_Rig = ((xRigU^a - xRigL^a)/(xRigU - xRigL));
        HRig(1,X.i) = -concv_slope_Rig;
        HRig(1,Z.i) = 1; 
        kRig(1,1) = -concv_slope_Rig*xRigL + xRigL^a;

        % Convex lowerbound
        % z > diff(x^a,x)|(x=xL)*(x - xL) + xL^a
        % z > diff(x^a,x)|(x=xmid)*(x - xmid) + xmid^a
        % z > diff(x^a,x)|(x=xU)*(x - xU) + xU^a
        HRig(2,X.i) = a*xRigL^(a-1);
        HRig(3,X.i) = a*xRigmid^(a-1);
        HRig(4,X.i) = a*xRigU^(a-1);
        HRig(2,Z.i) = -1;
        HRig(3,Z.i) = -1;
        HRig(4,Z.i) = -1;
        kRig(2,1) = a*xRigL^(a-1)*xRigL - xRigL^a;
        kRig(3,1) = a*xRigmid^(a-1)*xRigmid - xRigmid^a;
        kRig(4,1) = a*xRigU^(a-1)*xRigU - xRigU^a;      
        
        
        
        
        % Finding the vertices of the left concave bound and the right
        % convex bound
        % 4 vertices for each because we are using xL, xU and the midpoint
        LeftConcvVert = zeros(2,4);
        LeftConcvVert(:,1) = [xLefL; xLefL^a]; % x = xL, z = xL^a;
        LeftConcvVert(:,2) = HLef(1:2,[X.i,Z.i])\kLef(1:2,1);
        LeftConcvVert(:,3) = HLef(2:3,[X.i,Z.i])\kLef(2:3,1);
        LeftConcvVert(:,4) = [xLefU; xLefU^a]; % x = 0, z = 0^a;
        
        RightConvxVert = zeros(2,4);
        RightConvxVert(:,1) = [xRigL; xRigL^a]; % x = 0, z = 0^a;
        RightConvxVert(:,2) = HRig(2:3,[X.i,Z.i])\kRig(2:3,1);
        RightConvxVert(:,3) = HRig(3:4,[X.i,Z.i])\kRig(3:4,1);
        RightConvxVert(:,4) = [xRigU; xRigU^a]; % x = xU, z = xU^a;
        
                
        % Final envelope
        
        % Concave upperbound
        % z < ((xU^a - xL^a)/(xU - xL))*(x - xL) + xL^a        
        for j=4:-1:1
            concaveslope_candidate = (xU^a - LeftConcvVert(2,j))/(xU - LeftConcvVert(1,j));
            if(j==1) % There is no slope left to test, get the concave upperbound
                break;
            end            
            if(concaveslope_candidate < abs(HLef(j-1,X.i))) % Will give a concave upperbound
                break;
            end
        end
        concaveslope = concaveslope_candidate;
        Hconcv(1,X.i) = -concaveslope;
        Hconcv(1,Z.i) = 1; 
        kconcv(1,1) = -concaveslope_candidate*LeftConcvVert(1,j) + LeftConcvVert(2,j);               

        for k=1:j-1
            Hconcv = [Hconcv; HLef(k,:)];
            kconcv = [kconcv; kLef(k,:)];
        end        
        
        % Convex lowerbound
        % z > ((xU^a - xL^a)/(xU - xL))*(x - xL) + xL^a        
        for j=1:1:4
            convexslope_candidate = (RightConvxVert(2,j) - xL^a)/(RightConvxVert(1,j) - xL);
            if(j==4) % There is no slope left to compare with, get the convex lowerbound
                break;
            end
            if(convexslope_candidate < abs(HRig(j+1,X.i))) % Will give a convex lowerbound
                break;
            end
        end  
        convexslope = convexslope_candidate;
        Hconvx(1,X.i) = convexslope;
        Hconvx(1,Z.i) = -1; 
        kconvx(1,1) = convexslope_candidate*xL - xL^a;               

        for k=j+1:4
            Hconvx = [Hconvx; HRig(k,:)];
            kconvx = [kconvx; kRig(k,:)];
        end        
        
        
        % Build final equalities to add to Polyrelax Hrep
        Hadd = [Hconcv; Hconvx];
        kadd = [kconcv; kconvx];
  
        % Debugging section
%         keyboard;
%         xsample = xL:0.1:xU;
%         zsample = xsample.^a;
%         PLef = Polyhedron('A',HLef,'b',kLef);
%         PRig = Polyhedron('A',HRig,'b',kRig);
%         figure
%         plot(xsample,zsample,'LineWidth',2);
%         hold on
%         plot(PLef,'Color','b','Alpha',0.1);
%         plot(PRig,'Color','r','Alpha',0.1);  
%         for j=1:4
%             plot(LeftConcvVert(1,j),LeftConcvVert(2,j),'bx','MarkerSize',10);
%             plot(RightConvxVert(1,j),RightConvxVert(2,j),'rx','MarkerSize',10);
%         end
%         axis1 = axis;
%         Pconvx = Polyhedron('A',Hconvx,'b',kconvx,'LB',[axis1(1);axis1(3)],'UB',[axis1(2);axis1(4)]);
%         Pconcv = Polyhedron('A',Hconcv,'b',kconcv,'LB',[axis1(1);axis1(3)],'UB',[axis1(2);axis1(4)]);
%         plot(Pconvx,'Alpha',0);
%         plot(Pconcv,'Alpha',0);
%         axis(axis1);
%         figure;
%         plot(xsample,zsample,'LineWidth',2);
%         hold on
%         plot(Polyhedron('A',[Hconcv;Hconvx],'b',[kconcv;kconvx]),'Color','m','Alpha',0.1);        
%         keyboard;
        %error('For odd exponents, this function has not been implemented yet for arguments containing zero.'); 
        
        else
            
            % By finding the roots of the polynomials
            
            % xCV and xCC satisfy, respectively:
            %(n - 1)*(xCV)^n - n*xL(xCV)^(n-1) + (xL)^n = 0.
            %(n - 1)*(xCC)^n - n*xU(xCC)^(n-1) + (xU)^n = 0.
            
%             %xCV and xCC by roots function
%             xCV = roots([a-1,-a*xL,zeros(1,a-2),xL^a]);
%             xCC = roots([a-1,-a*xU,zeros(1,a-2),xU^a]);
%             
%             % Get only the positive real xCV and negative real xCC
%             xCV = xCV(imag(xCV)==0);
%             xCV = xCV(xCV>=0);
%             xCC = xCC(imag(xCC)==0);
%             xCC = xCC(xCC<=0);    

            % xCV and xCC by bisection method
            [xCV,xCC] = rootsbybisection(xL,xU,a);
            
            
            % Middle ("main") interval
            xLef = max(xL,xCC); % lower bound of main interval
            xRig = min(xU,xCV); % upper bound of main interval
           
            % Concave and convex secant slopes
            convx_slope = ((xRig^a - xL^a)/(xRig - xL));
            concv_slope = ((xU^a - xLef^a)/(xU - xLef));
            
            % Hrep of secants
            Hsec = zeros(2,Z.i); 
            ksec = zeros(2,1);
            
            % Convex secant
            % z > ((xU^a - xL^a)/(xU - xL))*(x - xL) + xL^a
            Hsec(1,X.i) = convx_slope;
            Hsec(1,Z.i) = -1; 
            ksec(1,1) = convx_slope*xL - xL^a;           
            
            % Concave secant
            % z < ((xU^a - xL^a)/(xU - xL))*(x - xL) + xL^a
            Hsec(2,X.i) = -concv_slope;
            Hsec(2,Z.i) = 1; 
            ksec(2,1) = -concv_slope*xLef + xLef^a;
            

            % Linear approximation for Right interval
            if(xU<=xCV) % If the Right interval does not exist
                HRig = zeros(0,Z.i);
                kRig = zeros(0,1);
            else % If it exists, linearize in 2 points (the first point would be equivalent to the convex secant)
                HRig = zeros(2,Z.i);
                kRig = zeros(2,1);               
                xRigM = 0.5*(xCV+xU);
                
                % Linearized convex lowerbound for Right interval
                % z > diff(x^a,x)|(x=xmid)*(x - xmid) + xmid^a
                % z > diff(x^a,x)|(x=xU)*(x - xU) + xU^a
                HRig(1,X.i) = a*xRigM^(a-1);
                HRig(2,X.i) = a*xU^(a-1);
                HRig(1,Z.i) = -1;
                HRig(2,Z.i) = -1;
                kRig(1,1) = a*xRigM^(a-1)*xRigM - xRigM^a;
                kRig(2,1) = a*xU^(a-1)*xU - xU^a;
            end
            
            % Linear approximation for Left interval
            if(xL>=xCC) % If the Left interval does not exist
                HLef = zeros(0,Z.i);
                kLef = zeros(0,1);
            else % If it exists, linearize in 2 points (the last, third point would be equivalent to the concave secant)
                HLef = zeros(2,Z.i);
                kLef = zeros(2,1);               
                xLefM = 0.5*(xL+xCC);
                
                % Linearized concave upper bound for Left interval
                % z < diff(x^a,x)|(x=xL)*(x - xL) + xL^a
                % z < diff(x^a,x)|(x=xmid)*(x - xmid) + xmid^a
                HLef(1,X.i) = -a*xL^(a-1);
                HLef(2,X.i) = -a*xLefM^(a-1);
                HLef(1,Z.i) = 1;
                HLef(2,Z.i) = 1;
                kLef(1,1) = -a*xL^(a-1)*xL + xL^a;
                kLef(2,1) = -a*xLefM^(a-1)*xLefM + xLefM^a;
            end    
            
            
            % Get the polyhedral enclosure Qj
            Hadd = [Hsec; HRig; HLef];
            kadd = [ksec; kRig; kLef];
                
            
        end
        
    end
    
        
    % Adds the inequality constraints to Polyrelax Hrep
    Polyrelax.addHrepHrow(Hadd,kadd);        
    
end

function power_noninteger(X,Z,a)

    xL = inf(X.x);
    xU = sup(X.x);
    xmid = 0.5*(xL+xU);
    Hadd = zeros(4,Z.i);
    kadd = zeros(4,1);
    
    if(xL<=0)
        error('For non-integer power, only positive arguments are currently implemented.');
    end
    
    if(a > 1) % Convex function for x > 0
  
        % Concave upperbound
        % z < ((xU^a - xL^a)/(xU - xL))*(x - xL) + xL^a
        concv_slope = ((xU^a - xL^a)/(xU - xL));
        Hadd(1,X.i) = -concv_slope;
        Hadd(1,Z.i) = 1; 
        kadd(1,1) = -concv_slope*xL + xL^a;

        % Convex lowerbound
        % z > diff(x^a,x)|(x=xL)*(x - xL) + xL^a
        % z > diff(x^a,x)|(x=xmid)*(x - xmid) + xmid^a
        % z > diff(x^a,x)|(x=xU)*(x - xU) + xU^a
        Hadd(2,X.i) = a*xL^(a-1);
        Hadd(3,X.i) = a*xmid^(a-1);
        Hadd(4,X.i) = a*xU^(a-1);
        Hadd(2,Z.i) = -1;
        Hadd(3,Z.i) = -1;
        Hadd(4,Z.i) = -1;
        kadd(2,1) = a*xL^(a-1)*xL - xL^a;
        kadd(3,1) = a*xmid^(a-1)*xmid - xmid^a;
        kadd(4,1) = a*xU^(a-1)*xU - xU^a;
        
    else % Concave function for x > 0
        
        % Convex lowerbound
        % z > ((xU^a - xL^a)/(xU - xL))*(x - xL) + xL^a
        concv_slope = ((xU^a - xL^a)/(xU - xL));
        Hadd(1,X.i) = concv_slope;
        Hadd(1,Z.i) = -1; 
        kadd(1,1) = concv_slope*xL - xL^a;

        % Concave upperbound
        % z < diff(x^a,x)|(x=xL)*(x - xL) + xL^a
        % z < diff(x^a,x)|(x=xmid)*(x - xmid) + xmid^a
        % z < diff(x^a,x)|(x=xU)*(x - xU) + xU^a
        Hadd(2,X.i) = -a*xL^(a-1);
        Hadd(3,X.i) = -a*xmid^(a-1);
        Hadd(4,X.i) = -a*xU^(a-1);
        Hadd(2,Z.i) = 1;
        Hadd(3,Z.i) = 1;
        Hadd(4,Z.i) = 1;
        kadd(2,1) = -a*xL^(a-1)*xL + xL^a;
        kadd(3,1) = -a*xmid^(a-1)*xmid + xmid^a;
        kadd(4,1) = -a*xU^(a-1)*xU + xU^a;
        
    end
    
    % Adds the inequality constraints to Polyrelax Hrep
    Polyrelax.addHrepHrow(Hadd,kadd);

end

function [xCV,xCC] = rootsbybisection(xL,xU,a)
% Find polynomial roots using fzero (which may use bisection)

persistent OPTIONS;
if(isempty(OPTIONS))
    %OPTIONS = optimset('Display','off');
    OPTIONS = optimset('fzero');
    OPTIONS.MyTOL = 1e-1;
end

bnds = max(abs([xL,xU]));

xCV = fzero(@(x) polconstraint(x,xL,a), [-OPTIONS.MyTOL,OPTIONS.MyTOL+bnds], OPTIONS);
xCC = fzero(@(x) polconstraint(x,xU,a), [-bnds-OPTIONS.MyTOL,OPTIONS.MyTOL], OPTIONS);

end

function out = polconstraint(x,xendp,a)
% Implements the polynomial constraint
out = (a-1)*x^a - a*xendp*x^(a-1) + xendp^a;
end

