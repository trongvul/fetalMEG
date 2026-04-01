% I adapted this code from Gaoxiang Ouyang (https://de.mathworks.com/matlabcentral/fileexchange/37289-permutation-entropy)

% Copyright (c) 2012, Gaoxiang
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
%
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


function [PermEn,PermEn_unnormed] = PermEn(y,m,t)
%  Calculate the permutation entropy
%  Input:   y: time series;
%           m: order of permuation entropy
%           t: delay time of permuation entropy, 
% Output: 
%           pe:    permuation entropy
%           hist:  the histogram for the order distribution
% Edited by Joel Frohlich 25.02.2022
%Ref: G Ouyang, J Li, X Liu, X Li, Dynamic Characteristics of Absence EEG Recordings with Multiscale Permutation %     %                             Entropy Analysis, Epilepsy Research, doi: 10.1016/j.eplepsyres.2012.11.003
%     X Li, G Ouyang, D Richards, Predictability analysis of absence seizures with permutation entropy, Epilepsy %     %                            Research,  Vol. 77pp. 70-74, 2007

warning('off', 'all')

ly = length(y);
permlist = perms(1:m);
c = nan(1,length(permlist));
    
 for j=1:ly-t*(m-1)
     [~,iv]=sort(y(j:t:j+t*(m-1)));
     for jj=1:length(permlist)
         if (abs(permlist(jj,:)-iv))==0
             if isnan(c(jj))
                 c(jj) = 1;
             else
                 c(jj) = c(jj) + 1;
             end
         end
     end
 end
%hist = c;
 
p = c/nansum(c);
pe = -nansum(p .* log(p));
PermEn = pe/log(factorial(m));
PermEn_unnormed = pe;

end
