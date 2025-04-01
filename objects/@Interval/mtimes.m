function ret = mtimes(a,b)
%* or MTIMES returns the product of two intervals

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

if(numel(a)>1 || numel(b)>1)
    % Matrix case
    flagtranspose = 0;
    if(~isnumeric(b))
        if(~isnumeric(a))
            if(size(a,2)~=size(b,1))
                error('Input dimensions mismatch.')
            else                      
                % If two interval matrices, resort to brute force multiplication to avoid conservatism. Can be slow.
                dim1 = size(a,1);
                dim2 = size(b,2);
                dimk = size(a,2);        
                ret = Interval.zeros(dim1,dim2);
                for i=1:dim1
                    for j=1:dim2
                        for k=1:dimk
                            ret(i,j) = ret(i,j) + a(i,k)*b(k,j);
                        end
                    end
                end 
            end
            return;
        else
            flagtranspose = 1;
        end
    end
    if(size(a,2)~=size(b,1))
        error('Input dimensions mismatch.')
    else
        % Faster implementation based on Algorithm 9.1 from Rump's Acta Automatica paper (2010)
        if(flagtranspose) % If the first argument is real and the second is interval, use c = (b.'*a.').'
            b_ = a.';
            a_ = b.';
            Ainf = inf(a_);
            Asup = sup(a_);
            mA = 0.5*(Ainf+Asup);
            rA = mA - Ainf;
            rC = rA*abs(b_);
            Csup = mA*b_ + rC;
            Cinf = mA*b_ - rC;
            ret = Interval(Cinf,Csup).';
        else % The first argument is interval and the second is real, use c = a*b
            Ainf = inf(a);
            Asup = sup(a);
            mA = 0.5*(Ainf+Asup);
            rA = mA - Ainf;
            rC = rA*abs(b);
            Csup = mA*b + rC;
            Cinf = mA*b - rC;
            ret = Interval(Cinf,Csup);           
        end
    end
                   
else
    % Scalar case
    a = Interval(a);
    b = Interval(b);

    S = [a.LB*b.LB, a.LB*b.UB, a.UB*b.LB, a.UB*b.UB];

    ret = Interval(min(S), max(S));
    
end

end