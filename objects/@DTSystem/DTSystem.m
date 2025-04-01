classdef DTSystem
% CLASS DTSYSTEM
% Defines a discrete-time system object, properties, and methods
%
% CONSTRUCTOR
%       Z = DTSystem(subclass,args); - Creates a discrete-time system of a given subclass
%
% INPUT
%   subclass: char array defining the subclass of the system (see the corresponding property)
%       args: for linear systems: pair-wise arguments defining matrices A, Bu, Bw, C, Du, Dv
%             for descriptor systems: pair-wise arguments defining matrices E, A, Bu, Bw, C, Du, Dv
%             for nonlinear systems: name (char array) of the system to be used in internal function handles
%
% SYNTAX EXAMPLES
%       DTSystem('linear','A',eye(3),'Bw',ones(3,2))
%       DTSystem('descriptor','A',eye(3),'Bw',ones(3,2),'E',[0 1 0; 0 0 1; 0 0 0])
%       DTSystem('nonlinear','Model1')
%       
% PROPERTIES
%       subclass: linear, descriptor, or nonlinear
%       E, A, Bu, Bw, C, Du, Dv: system matrices (for linear and descriptor systems)
%       f, g, h: function handles (for nonlinear systems)
%       nx, nu, ny, nw, nv: system dimensions (for nonlinear systems)
%       Ts: sampling time (for nonlinear systems, optional)
%
% METHODS
%       see help of each method

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
       
    properties (SetAccess=private)
        % Type: linear, descriptor, nonlinear
        subclass;
        % System matrices (for linear and descriptor systems)
        E; A; Bu; Bw; C; Dv;
        % Function handle structs (for nonlinear systems)
        f; % Dynamics 
        g; % Measurement
        h; % Invariants
        % System dimensions
        nx; nu; ny; nw; nv;        
        % Sampling time (for nonlinear systems, optional)
        Ts;
        % Last simulation data
        simdata;
        % Transformed space state variables (internal use, descriptor systems only)
        SVDvars;
    end
    
    methods

        function obj = DTSystem(subclass,varargin)
            if nargin <= 1 
                error('Not enough input arguments.');
            else
                obj.subclass = subclass;
                switch subclass
                    case 'linear'
                        
                        % Parse inputs
                        inp = inputParser;
                        inp.KeepUnmatched = false;
                        inp.addOptional('A',  [], @isnumeric);
                        inp.addOptional('Bu', [], @isnumeric);
                        inp.addOptional('Bw', [], @isnumeric);
                        inp.addOptional('C',  [], @isnumeric);
                        inp.addOptional('Dv', [], @isnumeric);
                        inp.parse(varargin{:});
                        inps = inp.Results;    
                        
                        [inps,dims] = DTSystem.checkinputs(inps,subclass);                        
                        
                        obj.A  = inps.A;
                        obj.Bu = inps.Bu;
                        obj.Bw = inps.Bw;
                        obj.C  = inps.C;
                        obj.Dv = inps.Dv;
                        obj.nx = dims.nx;
                        obj.nu = dims.nu;
                        obj.ny = dims.ny;
                        obj.nw = dims.nw;
                        obj.nv = dims.nv;                       
                        
                    case 'descriptor'
                        
                        % Parse inputs
                        inp = inputParser;
                        inp.KeepUnmatched = false;
                        inp.addOptional('E',  [], @isnumeric);                        
                        inp.addOptional('A',  [], @isnumeric);
                        inp.addOptional('Bu', [], @isnumeric);
                        inp.addOptional('Bw', [], @isnumeric);
                        inp.addOptional('C',  [], @isnumeric);
                        inp.addOptional('Dv', [], @isnumeric);
                        inp.parse(varargin{:});
                        inps = inp.Results;  
                        
                        [inps,dims] = DTSystem.checkinputs(inps,subclass);

                        obj.E  = inps.E;
                        obj.A  = inps.A;
                        obj.Bu = inps.Bu;
                        obj.Bw = inps.Bw;
                        obj.C  = inps.C;
                        obj.Dv = inps.Dv;
                        obj.nx = dims.nx;
                        obj.nu = dims.nu;
                        obj.ny = dims.ny;
                        obj.nw = dims.nw;
                        obj.nv = dims.nv;
                        
                        % Transformed model variables
                        [SVDvars.U,SVDvars.S,SVDvars.V] = svd(inps.E);
                        SVDvars.nof_nonzeros = sum(abs(diag(SVDvars.S))>0); % Number of nonzero eigenvalues
                        SVDvars.nof_zeros = size(inps.E,2) - SVDvars.nof_nonzeros;
                        SVDvars.Stilde = SVDvars.S(1:SVDvars.nof_nonzeros,1:SVDvars.nof_nonzeros);
                        A_SVD_ = (SVDvars.U\inps.A)/(SVDvars.V.');
                        B_SVD_ = SVDvars.U\inps.Bu;
                        Bw_SVD_ = SVDvars.U\inps.Bw;
                        SVDvars.Atilde = SVDvars.Stilde\A_SVD_(1:SVDvars.nof_nonzeros,:);
                        SVDvars.Btilde = SVDvars.Stilde\B_SVD_(1:SVDvars.nof_nonzeros,:);
                        SVDvars.Bwtilde = SVDvars.Stilde\Bw_SVD_(1:SVDvars.nof_nonzeros,:);
                        SVDvars.Ahat = A_SVD_(SVDvars.nof_nonzeros+1:end,:);
                        SVDvars.Bhat = B_SVD_(SVDvars.nof_nonzeros+1:end,:);
                        SVDvars.Bwhat = Bw_SVD_(SVDvars.nof_nonzeros+1:end,:);
                        SVDvars.T = inv(SVDvars.V.'); % Transformation matrix from the transformed space to the original space 
                        obj.SVDvars = SVDvars;
                        
                    case 'nonlinear'
                        name = varargin{1};
                        if(ischar(name))
                            obj.f = util.fhandle([name,'_f']);
                            obj.g = util.fhandle([name,'_g']);
                            obj.h = util.fhandle([name,'_h']);
                        else
                            error('Input name for nonlinear system constructor must be a char array.');
                        end
                        
                        % Check dimensions
                        nx = varargin{2};
                        nu = varargin{3};
                        ny = varargin{4};
                        nw = varargin{5};
                        nv = varargin{6};
                        if((~isnumeric(nx)||isempty(nx)||numel(nx)>1)||...
                           (~isnumeric(nu)||isempty(nu)||numel(nu)>1)||...
                           (~isnumeric(ny)||isempty(ny)||numel(ny)>1)||...
                           (~isnumeric(nw)||isempty(nw)||numel(nw)>1)||...
                           (~isnumeric(nv)||isempty(nv)||numel(nv)>1))
                            error('Input dimensions must be scalars.');
                        elseif(~((nx>=0)&&(nu>=0)&&(ny>=0)&&(nw>=0)&&(nv>=0)))
                            error('Input dimensions must be positive or zero.');
                        end
                        obj.nx = nx;
                        obj.nu = nu;
                        obj.ny = ny;
                        obj.nw = nw;
                        obj.nv = nv;                        
                        
                        % Sampling time
                        if(nargin>7)
                            Ts = varargin{7};
                            if(isnumeric(Ts)&&(numel(Ts)==1))
                                if(Ts > 0)
                                    obj.Ts = Ts;
                                else
                                    error('Input sampling time for nonlinear system constructor must be a positive scalar.');
                                end
                            else
                                error('Input sampling time for nonlinear system constructor must be a positive scalar.');
                            end
                        else
                            obj.Ts = [];
                        end
                    otherwise
                        error('Invalid subclass in DTSystem constructor.');
                end
            end
        end 
      
    end
    
    methods(Static, Access=private) 
        
        function [val,dims] = checkinputs(inps,subclass)
            % Input check for linear and descriptor systems
                        
            % Flags
            hasA  = ~isempty(inps.A);
            hasBu = ~isempty(inps.Bu);
            hasBw = ~isempty(inps.Bw);   
            hasC  = ~isempty(inps.C);
            hasDv = ~isempty(inps.Dv);  
                       
            % Number of system inputs
            if(hasBu)
                nu = size(inps.Bu,2);
            else
                nu = 0;
            end
                       
            % Dynamics check
            % Requires at least one of A, Bu, or Bw to be defined
            
            if((~hasA)&&(~hasBu)&&(~hasBw))
                error('For linear or descriptor systems, at least one of A, Bu, or Bw must be defined.');
            end
            
            if(hasA)
                nx = size(inps.A,1); % Number of states
                if(size(inps.A,2)~=nx)
                    error('Matrix A must be square.');
                end
                if(hasBu)
                    if(size(inps.Bu,1)~=nx)
                        DTSystem.errordimensions;
                    end
                else
                    inps.Bu = zeros(nx,nu);
                end
                if(hasBw)
                    if(size(inps.Bw,1)~=nx)
                        DTSystem.errordimensions;
                    end
                else
                    inps.Bw = zeros(nx,0);
                end  
            elseif(hasBu)
                nx = size(inps.Bu,1); % Number of states
                inps.A = zeros(nx,nx);
                if(hasBw)
                    if(size(inps.Bw,1)~=nx)
                        DTSystem.errordimensions;
                    end
                else
                    inps.Bw = zeros(nx,0);
                end  
            elseif(hasBw)
                nx = size(inps.Bw,1); % Number of states
                inps.A = zeros(nx,nx);
                inps.Bu = zeros(nx,nu);
            end
            nw = size(inps.Bw,2);
            
            if(strcmp(subclass,'descriptor'))
                hasE = ~isempty(inps.E);
                if(hasE)
                    if(((size(inps.E,1)~=size(inps.E,2))&&(size(inps.E,1)~=nx)))
                        DTSystem.errordimensions;
                    end
                else   
                    inps.E = eye(nx);
                end
            end
                       
            % Measurement check            
            % Will allow to not have any output matrix
            
            if(hasC)
                ny = size(inps.C,1); % Number of measurements
                if(size(inps.C,2)~=nx)
                    error('Matrix C does not match the number of states.');
                end
                if(hasDv)
                    if(size(inps.Dv,1)~=ny)
                        DTSystem.errordimensions;
                    end
                else
                    inps.Dv = zeros(ny,0);
                end  
            elseif(hasDv)
                ny = size(inps.Dv,1); % Number of measurements
                inps.C = zeros(ny,nx);
                inps.Du = zeros(ny,nu);
            else
                ny = 0;
                inps.C = zeros(0,nx);
                inps.Dv = [];
            end
            nv = size(inps.Dv,2);
            
            if(strcmp(subclass,'descriptor'))
                hasE = ~isempty(inps.E);
                if(hasE)
                    if(((size(inps.E,1)~=size(inps.E,2))&&(size(inps.E,1)~=nx)))
                        DTSystem.errordimensions;
                    end
                else   
                    inps.E = eye(nx);
                end
            end            
                       
            val = inps;
            dims = struct('nx',nx,'nu',nu,'ny',ny,'nw',nw,'nv',nv);
                
        end
        
        function errordimensions
            error('System matrices have incompatible dimensions.');
        end
        
    end
end