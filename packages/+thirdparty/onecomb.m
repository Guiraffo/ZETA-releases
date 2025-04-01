function comb = onecomb(N,K,M)
% ONECOMB - obtain the M-th combination from NCHOOSEK(1:N,K)
%   comb = ONECOMB(N,K,M) returns the m-th combination of the sorted 
%   list of all combinations from NCHOOSEK. 
%   N, K, M are non-negative scalar, comb has size 1-by-K.
%
%   See also PERMS, NCHOOSEK
%        and COMBN and ALLCOMB on the File Exchange
%
% Algorithm: For given N, K and M, where 1 <= M <= C(N,K), 
%            M2 = C(N,K)-M can be represented uniquely in the form
%            
%            M2 = C(m1,K) + C(m2,K-1) + ... + C(mK,1)
%            
%            where N > m1 > m2 > ... > mK. The ordered vector
%            (m1,m2,...,mK) is called the combinadic of M2 given N and K.
%            Within this scheme C(N,K) = 0 if N < K.
%            The M-th combination of NCHOOSEK is obtained from the
%            combinadic of M2. More information in the article cited 
%            below and references within.
% Reference: article by James McCaffrey (MSDN) accessed Feb09
% http://msdn.microsoft.com/en-us/library/aa289166(VS.71).aspx
% for Matlab (should work for most versions)
% version 1.0 (Feb 2009)
% (c) Darren Rowland
% email: darrenjrowland@hotmail.com
%
% Keywords: combinations, combinadic
error(nargchk(3,3,nargin)) ;
if numel(N) ~= 1 || N <= 0 || N ~= round(N)
  error('onecomb:InvalidArg1',...
        'The first input has to be a non-negative integer');
end
if numel(M) ~= 1 || M <= 0 || M ~= round(M)
  error('onecomb:InvalidArg3',...
        'The third input has to be a non-negative integer');
end
NCK = nchoosek(N,K);
if M>NCK 
     error('onecomb:largeM','M should not exceed N!/(K!(N-K)!)');
end
comb = zeros(1,K);
K2 = K;
N2 = N;
M2 = NCK - M;
for ii=1:K
    if M2 > 0
        % Find largest N such that C(N,K) <= M2
        while(NCK>M2)
            % Reducing N by 1, therefore require C(N-1,K)
            % C(N-1,K) = C(N,K)/N*(N-K)
            NCK = NCK/N2*(N2-K2);
            N2 = N2 - 1;
        end
        M2 = M2 - NCK;
        comb(ii) = N2;
        if ii == K
            % combinadic is complete, short circuit to avoid rare
            % divide-by-zero in the next line
            break
        end
        % Reducing K and N by 1, therefore require C(N-1,K-1)
        % C(N-1,K-1) = C(N,K)/N*K
        NCK = NCK/N2*K2;
        N2 = N2 - 1;
        K2 = K2 - 1;
    else
        % When M2 reaches zero, the remaining values in comb can 
        % be filled in automatically
        comb(ii:K) = K2-1:-1:0;
        break
    end
end
% Convert combinadic into valid combination
comb = N - comb;