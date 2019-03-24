%--------------------------------------------------------------------------
% function to try to filter the output of the wing averaging fits. performs
% a hampel filter for outlier removal and then a low pass filter for
% smoothing
%--------------------------------------------------------------------------
function vecs_filt = filter_wing_vecs(vecs, hampel_k, hampel_sigma, ...
    filt_lvl,normFlag) 

if ~exist('normFlag','var')
    normFlag = false ;
end

% filter design 
d1 = designfilt('lowpassiir','FilterOrder', 8, 'SampleRate', 8000,...
    'HalfPowerFrequency', filt_lvl ,'DesignMethod', 'butter') ;
% make sure the data is in a Nx3 matrix
if size(vecs,2) > size(vecs,1)
    vecs = vecs' ; 
end

% hampel filter
[~, hampel_idx] = hampel(vecs, hampel_k, hampel_sigma) ; 
vecs(hampel_idx) = nan ; 

%vecs_filt = vecs ; 

% interpolate nan values and then filter results
nan_idx = isnan(vecs) ; 
frames = (1:size(vecs,1))' ;
vecs_filt = nan(size(vecs)) ; 
for i = 1:size(vecs,2)
    c_fit = fit(frames(~nan_idx(:,i)), vecs(~nan_idx(:,i),i), 'cubicinterp') ;
    idx1 = find(~nan_idx(:,i),1,'first') ; 
    idx2 = find(~nan_idx(:,i),1,'last') ; 
    vecs_interp = c_fit(frames(idx1:idx2)) ; 
    vecs_filt(idx1:idx2,i) = filtfilt(d1,vecs_interp) ; 
end

if normFlag
   vecs_filt = vecs_filt./repmat(myNorm(vecs_filt),1,3) ;  
end
end