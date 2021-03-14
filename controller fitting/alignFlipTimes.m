% -------------------------------------------------------------------------
% solve assignment problem for wing flip times
% -------------------------------------------------------------------------
function [timesR_out, timesL_out, R_idx, L_idx] = ...
    alignFlipTimes(timesR, timesL, noMatchCost) 
% if not given, provide cost for not matching (note, here we're using
% seconds)
if ~exist('noMatchCost','var') || isempty(noMatchCost)
    if (mean(diff(timesR)) > 1)
        noMatchCost = 1 ;
    else
        noMatchCost = 1e-3 ; 
    end
end

% make sure the time entries are in column form 
if (size(timesR,1) < size(timesR,2))
    timesR = timesR' ; 
    timesL = timesL' ; 
end

% calculate distances between all cut off times and solve matching problem
distMat = pdist2(timesR, timesL) ; 
M = matchpairs(distMat, noMatchCost) ;
R_idx = M(:,1) ; 
L_idx = M(:,2) ; 
timesR_out = timesR(R_idx) ; 
timesL_out = timesL(L_idx) ; 

end