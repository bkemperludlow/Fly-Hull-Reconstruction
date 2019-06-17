%--------------------------------------------------------------------------
% function to pull out indices of when a signal goes above a threshold. the
% chunks of continuous indices are split up
%--------------------------------------------------------------------------
function idx_list = idx_by_thresh(signal, thresh)

if nargin < 2 
    thresh = 0.1 ; 
end

if size(signal,1) < size(signal,2)
    signal = signal' ; 
end

idxs = find(signal > thresh) ; 

if ~isempty(idxs)
    split_idxs = find(diff(idxs) > 1) ; 
else
    idx_list = [] ; 
    return ;
end

split_idxs = [0 ; split_idxs ; length(idxs)] ; 
idx_list = cell(length(split_idxs)-1,1) ; 

for i = 1:(length(split_idxs) - 1)
    idx_list{i} = idxs((split_idxs(i)+1):split_idxs(i+1)) ; 
end

not_empty_ind = cellfun(@(x) ~isempty(x), idx_list) ;  
idx_list = idx_list(not_empty_ind) ; 

end