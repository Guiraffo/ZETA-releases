function disp(obj)
% displays an Interval

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

    row = size(obj,1);
    column = size(obj,2);
    
    %disp(' ');
    disp(['Interval with endpoints:'])
    %disp(' ');
    
    for i = 1:row   
        row_string= '     ';
        
        for j = 1:column
            %row_string = [row_string, '[ ', num2str(obj(i,j).LowerBound), ' , ', num2str(obj(i,j).UpperBound), ' ]', '    ']; 
            row_string = [row_string, '[ ', num2str(obj(i,j).LB), ', ', num2str(obj(i,j).UB), ']', '  '];          
        end
        disp(row_string);
    end
    %disp(' ');
    
end
