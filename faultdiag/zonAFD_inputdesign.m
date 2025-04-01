function [utilde, deltahat_m, cputimes] = zonAFD_inputdesign(sysarr,r,s,X0,W,V,U,OPTIONS)
%ZONAFD_INPUTDESIGN design a sequence of separating inputs for a given 
%                   array of linear models and separation horizon
%
%   SYNTAX: utilde = ZONAFD_INPUTDESIGN(A,B,Dw,r,C,Dv,s,X0,W,V,N)
%
%   INPUTS
%      sysarr: array of linear models (DTSystem objects)
%           r: cell array of process bias
%           s: cell array of measurement bias
%          X0: initial zonotope X0
%           W: uncertainty zonotope W
%           V: uncertainty zonotope V
%           U: set of admissible inputs (struct with fields H and k 
%              describing a convex polyhedron in halfspace representation)
%     OPTIONS: options structure built using util.faultdiagopts (see help)
%              
%   OUTPUT
%      utilde: separating input sequence
%              

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

% Get config variables from struct
R = OPTIONS.R;
N = OPTIONS.N;
max_order = OPTIONS.max_order;
epsilon = OPTIONS.epsilon;

A  = {sysarr.A};
Bu = {sysarr.Bu};
Bw = {sysarr.Bw};
C  = {sysarr.C};
Dv = {sysarr.Dv};

nof_models = length(sysarr);
nof_inps = size(Bu{1},2);

% Get admissible set for the input sequence
[Utilde.H, Utilde.k] = util.repmatpoly(U.H,U.k,N);

% Get augmented variables

Psi = cell(nof_models,1);
Tildevars = cell(nof_models,1);
for j=1:nof_models
   [~,Psi{j},Tildevars{j}] = AFD.LTI_reachable_sets(r{j},A{j},Bu{j},Bw{j},s{j},C{j},Dv{j},X0,W,V,zeros(nof_inps,N));
end


% Generate separating constraints

Ntildeq = cell(nof_models,nof_models);
Ztildeq = cell(nof_models,nof_models);
Gq = cell(nof_models,nof_models);
cq = cell(nof_models,nof_models);
for i=1:nof_models-1
    for j=i+1:nof_models
        
        % Pair combination variables
        Ntildeq{i,j} = Tildevars{j}.C*Tildevars{j}.Bu - Tildevars{i}.C*Tildevars{i}.Bu;
        Ztildeq{i,j} = Zonotope(Psi{i}.c - Psi{j}.c, [Psi{i}.G, Psi{j}.G]);

        % Order reduction
        if(~isempty(max_order))
           if((size(Ztildeq{i,j}.G,2)/size(Ztildeq{i,j}.G,1)) > max_order)                
               ngmax = max_order*size(Ztildeq{i,j}.G,1);
               Ztildeq{i,j} = reduction(Ztildeq{i,j},ngmax);
           end       
        end        
       
        Gq{i,j} = Ztildeq{i,j}.G;
        cq{i,j} = Ztildeq{i,j}.c;             
        
    end
end



% OPTIONS
OPTIONS.Display = 'off';

space_dim = size(Gq{1,2},1);
nof_gens = size(Gq{1,2},2);

Rtildeaux = repmat({R},N,1);
Rtilde = blkdiag(Rtildeaux{:});

% Length of the input sequence
usize = nof_inps*N;

% Number of halfspaces in Utilde
nof_h_utilde = size(Utilde.H,1);

% Number of possible pair-wise combinations of models
nof_q = nchoosek(nof_models,2);

% Find the upper bounds deltahat_m
tic;
deltahat_m = AFD.LTI_find_upperbounds(Utilde,Ntildeq,Ztildeq,N,OPTIONS.deltahat_guess,[]);
cputimes.findbounds = toc;


% Builds and solves the mixed-integer program

switch OPTIONS.functional
    case 'quadratic'
    
% Mixed-integer quadratic program

Hess = blkdiag(2*Rtilde, zeros(nof_q*(2 + 5*nof_gens + space_dim)));
f = [                              zeros(usize,1);
     zeros(nof_q*(2 + 5*nof_gens + space_dim),1)];
Aineq = [  Utilde.H, zeros(nof_h_utilde,nof_q*(2 + 5*nof_gens + space_dim))];
bineq =    Utilde.k;
Aeq = [];
beq = [];    
ctype = repmat('C',1,usize);

q = 1;
for i=1:nof_models-1
    for j=i+1:nof_models

        Aineq = [Aineq; 
                 zeros(1, usize),                zeros(1, (q-1)*(2 + 5*nof_gens + space_dim)),                -1, 0,  zeros(1, nof_gens), zeros(1, space_dim), zeros(1, nof_gens), zeros(1, nof_gens), zeros(1, nof_gens), zeros(1, nof_gens),                                                                                                                    zeros(1, (nof_q-q)*(2 + 5*nof_gens + space_dim)); 
                 zeros(1, usize),                zeros(1, (q-1)*(2 + 5*nof_gens + space_dim)),                 1, 0,  zeros(1, nof_gens), zeros(1, space_dim), zeros(1, nof_gens), zeros(1, nof_gens), zeros(1, nof_gens), zeros(1, nof_gens),                                                                                                                    zeros(1, (nof_q-q)*(2 + 5*nof_gens + space_dim));             
                 zeros(1, usize),                zeros(1, (q-1)*(2 + 5*nof_gens + space_dim)),                -1, 1,  zeros(1, nof_gens), zeros(1, space_dim), zeros(1, nof_gens), zeros(1, nof_gens), zeros(1, nof_gens), zeros(1, nof_gens),                                                                                                                    zeros(1, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    zeros(nof_gens,1),  -ones(nof_gens,1),      eye(nof_gens), zeros(nof_gens, space_dim),  zeros(nof_gens, nof_gens), zeros(nof_gens, nof_gens), zeros(nof_gens, nof_gens), zeros(nof_gens, nof_gens),      zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    zeros(nof_gens,1),  -ones(nof_gens,1),     -eye(nof_gens), zeros(nof_gens, space_dim),  zeros(nof_gens, nof_gens), zeros(nof_gens, nof_gens), zeros(nof_gens, nof_gens), zeros(nof_gens, nof_gens),      zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    zeros(nof_gens,1),  zeros(nof_gens,1),    zeros(nof_gens), zeros(nof_gens, space_dim),  -eye(nof_gens), zeros(nof_gens), zeros(nof_gens), zeros(nof_gens),                                                                       zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    zeros(nof_gens,1),  zeros(nof_gens,1),    zeros(nof_gens), zeros(nof_gens, space_dim), zeros(nof_gens),  -eye(nof_gens), zeros(nof_gens), zeros(nof_gens),                                                                       zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    zeros(nof_gens,1),  zeros(nof_gens,1),    zeros(nof_gens), zeros(nof_gens, space_dim),   eye(nof_gens), zeros(nof_gens),  -eye(nof_gens), zeros(nof_gens),                                                                       zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    zeros(nof_gens,1),  zeros(nof_gens,1),    zeros(nof_gens), zeros(nof_gens, space_dim), zeros(nof_gens),   eye(nof_gens), zeros(nof_gens),  -eye(nof_gens),                                                                       zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),     ones(nof_gens,1),  zeros(nof_gens,1),     -eye(nof_gens), zeros(nof_gens, space_dim), zeros(nof_gens), zeros(nof_gens), 2*(1+deltahat_m{i,j})*eye(nof_gens),                     zeros(nof_gens),                                zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),     ones(nof_gens,1),  zeros(nof_gens,1),      eye(nof_gens), zeros(nof_gens, space_dim), zeros(nof_gens), zeros(nof_gens),                     zeros(nof_gens), 2*(1+deltahat_m{i,j})*eye(nof_gens),                                zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    -ones(nof_gens,1),  zeros(nof_gens,1),      eye(nof_gens), zeros(nof_gens, space_dim), zeros(nof_gens), zeros(nof_gens),  zeros(nof_gens), zeros(nof_gens),                                                                       zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    -ones(nof_gens,1),  zeros(nof_gens,1),     -eye(nof_gens), zeros(nof_gens, space_dim), zeros(nof_gens), zeros(nof_gens),  zeros(nof_gens), zeros(nof_gens),                                                                       zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim))];

        bineq = [bineq;
                            -epsilon;
                     deltahat_m{i,j};            
                                   1;
                   zeros(nof_gens,1);
                   zeros(nof_gens,1);
                   zeros(nof_gens,1);
                   zeros(nof_gens,1);
                   zeros(nof_gens,1);
                   zeros(nof_gens,1);      
             (-1 + 2*(1+deltahat_m{i,j}))*ones(nof_gens,1);
             (-1 + 2*(1+deltahat_m{i,j}))*ones(nof_gens,1);
                    ones(nof_gens,1);
                    ones(nof_gens,1)];            


        Aeq = [Aeq;
                         Ntildeq{i,j}, zeros(space_dim, (q-1)*(2 + 5*nof_gens + space_dim)), zeros(space_dim,1), zeros(space_dim,1),            -Gq{i,j},    zeros(space_dim), zeros(space_dim,nof_gens), zeros(space_dim,nof_gens), zeros(space_dim,nof_gens), zeros(space_dim,nof_gens),  zeros(space_dim, (nof_q-q)*(2 + 5*nof_gens + space_dim));
               zeros(nof_gens, usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),  zeros(nof_gens,1),  zeros(nof_gens,1),     zeros(nof_gens),           Gq{i,j}.',            -eye(nof_gens),             eye(nof_gens),           zeros(nof_gens),           zeros(nof_gens),   zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                      zeros(1, usize),         zeros(1, (q-1)*(2 + 5*nof_gens + space_dim)),                  0,                  0,  zeros(1, nof_gens), zeros(1, space_dim),         ones(1, nof_gens),         ones(1, nof_gens),         zeros(1,nof_gens),         zeros(1,nof_gens),          zeros(1, (nof_q-q)*(2 + 5*nof_gens + space_dim))];             


        beq = [beq;
                           cq{i,j};
                 zeros(nof_gens,1);
                                 1];   

        ctype = strcat([ctype,repmat('C',1,1+1+nof_gens+space_dim+2*nof_gens), repmat('B',1,2*nof_gens)]); 
        q = q + 1;

    end
end

tic
[x,~,exitflag] = optim.solvemiqp(Hess,f,Aineq,bineq,Aeq,beq,ctype,OPTIONS);
cputimes.inputdesign = toc;

if(ischar(exitflag)) % Gurobi
    if(strcmp(exitflag,'OPTIMAL')) % Success
        utilde = x(1:usize);
    else % Not success 
        error(strcat(['Exit flag ',exitflag,' in zonAFD_inputdesign']));
    end    
elseif(exitflag==1) % intlinprog: Success
    utilde = x(1:usize);  
else % Not success 
    error(strcat(['Exit flag ',num2str(exitflag),' in zonAFD_inputdesign.']));
end    

    case 'linear'

% Mixed-integer linear program (1-norm)
        
tsize = usize; % size of the auxiliary decision variable t (1-norm)        

f = [   ones(tsize,1);
       zeros(usize,1);
     zeros(nof_q*(2 + 5*nof_gens + space_dim),1)];
Aineq = [  zeros(nof_h_utilde, tsize), Utilde.H, zeros(nof_h_utilde,nof_q*(2 + 5*nof_gens + space_dim))];
bineq =    Utilde.k;
Aeq = [];
beq = [];    
ctype = repmat('C',1,tsize + usize);


% 1-norm constraints
Aineq = [Aineq; 
         -eye(tsize),  Rtilde, zeros(tsize, nof_q*(2 + 5*nof_gens + space_dim));
         -eye(tsize), -Rtilde, zeros(tsize, nof_q*(2 + 5*nof_gens + space_dim))];
     
bineq = [bineq; 
         zeros(tsize,1);
         zeros(tsize,1)];
     
q = 1;
for i=1:nof_models-1
    for j=i+1:nof_models

        Aineq = [Aineq; 
                 zeros(1, tsize+usize),                zeros(1, (q-1)*(2 + 5*nof_gens + space_dim)),                -1, 0,  zeros(1, nof_gens), zeros(1, space_dim), zeros(1, nof_gens), zeros(1, nof_gens), zeros(1, nof_gens), zeros(1, nof_gens),                                                                                                                    zeros(1, (nof_q-q)*(2 + 5*nof_gens + space_dim)); 
                 zeros(1, tsize+usize),                zeros(1, (q-1)*(2 + 5*nof_gens + space_dim)),                 1, 0,  zeros(1, nof_gens), zeros(1, space_dim), zeros(1, nof_gens), zeros(1, nof_gens), zeros(1, nof_gens), zeros(1, nof_gens),                                                                                                                    zeros(1, (nof_q-q)*(2 + 5*nof_gens + space_dim));             
                 zeros(1, tsize+usize),                zeros(1, (q-1)*(2 + 5*nof_gens + space_dim)),                -1, 1,  zeros(1, nof_gens), zeros(1, space_dim), zeros(1, nof_gens), zeros(1, nof_gens), zeros(1, nof_gens), zeros(1, nof_gens),                                                                                                                    zeros(1, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, tsize+usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    zeros(nof_gens,1),  -ones(nof_gens,1),      eye(nof_gens), zeros(nof_gens, space_dim),  zeros(nof_gens, nof_gens), zeros(nof_gens, nof_gens), zeros(nof_gens, nof_gens), zeros(nof_gens, nof_gens),      zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, tsize+usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    zeros(nof_gens,1),  -ones(nof_gens,1),     -eye(nof_gens), zeros(nof_gens, space_dim),  zeros(nof_gens, nof_gens), zeros(nof_gens, nof_gens), zeros(nof_gens, nof_gens), zeros(nof_gens, nof_gens),      zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, tsize+usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    zeros(nof_gens,1),  zeros(nof_gens,1),    zeros(nof_gens), zeros(nof_gens, space_dim),  -eye(nof_gens), zeros(nof_gens), zeros(nof_gens), zeros(nof_gens),                                                                       zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, tsize+usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    zeros(nof_gens,1),  zeros(nof_gens,1),    zeros(nof_gens), zeros(nof_gens, space_dim), zeros(nof_gens),  -eye(nof_gens), zeros(nof_gens), zeros(nof_gens),                                                                       zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, tsize+usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    zeros(nof_gens,1),  zeros(nof_gens,1),    zeros(nof_gens), zeros(nof_gens, space_dim),   eye(nof_gens), zeros(nof_gens),  -eye(nof_gens), zeros(nof_gens),                                                                       zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, tsize+usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    zeros(nof_gens,1),  zeros(nof_gens,1),    zeros(nof_gens), zeros(nof_gens, space_dim), zeros(nof_gens),   eye(nof_gens), zeros(nof_gens),  -eye(nof_gens),                                                                       zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, tsize+usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),     ones(nof_gens,1),  zeros(nof_gens,1),     -eye(nof_gens), zeros(nof_gens, space_dim), zeros(nof_gens), zeros(nof_gens), 2*(1+deltahat_m{i,j})*eye(nof_gens),                     zeros(nof_gens),                                zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, tsize+usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),     ones(nof_gens,1),  zeros(nof_gens,1),      eye(nof_gens), zeros(nof_gens, space_dim), zeros(nof_gens), zeros(nof_gens),                     zeros(nof_gens), 2*(1+deltahat_m{i,j})*eye(nof_gens),                                zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, tsize+usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    -ones(nof_gens,1),  zeros(nof_gens,1),      eye(nof_gens), zeros(nof_gens, space_dim), zeros(nof_gens), zeros(nof_gens),  zeros(nof_gens), zeros(nof_gens),                                                                       zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                 zeros(nof_gens, tsize+usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),    -ones(nof_gens,1),  zeros(nof_gens,1),     -eye(nof_gens), zeros(nof_gens, space_dim), zeros(nof_gens), zeros(nof_gens),  zeros(nof_gens), zeros(nof_gens),                                                                       zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim))];

        bineq = [bineq;
                            -epsilon;
                     deltahat_m{i,j};            
                                   1;
                   zeros(nof_gens,1);
                   zeros(nof_gens,1);
                   zeros(nof_gens,1);
                   zeros(nof_gens,1);
                   zeros(nof_gens,1);
                   zeros(nof_gens,1);      
             (-1 + 2*(1+deltahat_m{i,j}))*ones(nof_gens,1);
             (-1 + 2*(1+deltahat_m{i,j}))*ones(nof_gens,1);
                    ones(nof_gens,1);
                    ones(nof_gens,1)];            


        Aeq = [Aeq;
               zeros(space_dim, tsize),            Ntildeq{i,j}, zeros(space_dim, (q-1)*(2 + 5*nof_gens + space_dim)), zeros(space_dim,1), zeros(space_dim,1),            -Gq{i,j},    zeros(space_dim), zeros(space_dim,nof_gens), zeros(space_dim,nof_gens), zeros(space_dim,nof_gens), zeros(space_dim,nof_gens),  zeros(space_dim, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                zeros(nof_gens, tsize),  zeros(nof_gens, usize),  zeros(nof_gens, (q-1)*(2 + 5*nof_gens + space_dim)),  zeros(nof_gens,1),  zeros(nof_gens,1),     zeros(nof_gens),           Gq{i,j}.',            -eye(nof_gens),             eye(nof_gens),           zeros(nof_gens),           zeros(nof_gens),   zeros(nof_gens, (nof_q-q)*(2 + 5*nof_gens + space_dim));
                       zeros(1, tsize),         zeros(1, usize),         zeros(1, (q-1)*(2 + 5*nof_gens + space_dim)),                  0,                  0,  zeros(1, nof_gens), zeros(1, space_dim),         ones(1, nof_gens),         ones(1, nof_gens),         zeros(1,nof_gens),         zeros(1,nof_gens),          zeros(1, (nof_q-q)*(2 + 5*nof_gens + space_dim))];             


        beq = [beq;
                           cq{i,j};
                 zeros(nof_gens,1);
                                 1];   

        ctype = strcat([ctype,repmat('C',1,1+1+nof_gens+space_dim+2*nof_gens), repmat('B',1,2*nof_gens)]); 
        q = q + 1;

    end
end

tic
[x,~,exitflag] = optim.solvemilp(f,Aineq,bineq,Aeq,beq,ctype,OPTIONS);
cputimes.inputdesign = toc;
     
 
if(ischar(exitflag)) % Gurobi
    if(strcmp(exitflag,'OPTIMAL')) % Success
        utilde = x(tsize+1:tsize+usize);
    else % Not success 
        error(strcat(['Exit flag ',exitflag,' in zonAFD_inputdesign']));
    end    
elseif(exitflag==1) % intlinprog: Success
    utilde = x(tsize+1:tsize+usize);  
else % Not success 
    error(strcat(['Exit flag ',num2str(exitflag),' in zonAFD_inputdesign.']));
end    

    otherwise
        error('Invalid choice of cost function for input design.');
end
        




end

