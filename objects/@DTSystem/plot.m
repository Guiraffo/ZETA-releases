function plot(sys)
%PLOT the last simulation results of a discrete-time system
%
%   SYNTAX: PLOT(sys)
%
%   INPUTS
%         sys: discrete-time system (dtsystem)

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

simdata = sys.simdata;
if(isempty(simdata))
    error('No simulation data found!');
end

nx = sys.nx;
ny = sys.ny;
nu = sys.nu;
Ts = sys.Ts;
simlength = size(simdata.x,2);

if(isempty(Ts))
    simtime = 0:1:(simlength-1);
    timelabel = 'k';
else
    simtime = 0:Ts:(simlength-1)*Ts;
    timelabel = 'time';
end
    

% Plot states
for j=1:nx
    figure;
    hold on;
    plot(simtime,simdata.x(j,:),'LineWidth',1);
    legend(['x_',num2str(j)]);
    xlabel(timelabel);
    ylabel('x');
    box on;
end

% Plot measurements, if any
for j=1:ny
    figure;
    hold on;
    plot(simtime,simdata.y(j,:),'r','LineWidth',1);
    legend(['y_',num2str(j)]);
    xlabel(timelabel);
    ylabel('y');
    box on;
end

% Plot known inputs, if any
for j=1:nu
    figure;
    hold on;
    plot(simtime,simdata.u(j,:),'m','LineWidth',1);
    legend(['u_',num2str(j)]);
    xlabel(timelabel);
    ylabel('u');
    box on;
end

end
