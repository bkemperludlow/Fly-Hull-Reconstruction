% -------------------------------------------------------------------------
% function to use the vector between left and right wing centers of mass to
% estimate roll, find frames with left/right wing swap, and identify other
% poorly-reconstructed frames (hopefully)
%
% TO DO:
%   - make sure of transformations to/from body frame
%   - use smoothed roll to identify areas where we should swap the wings
%   - figure out best form of output for this function (data stucture?
%   rotation matrices? just roll? assign new roll vectors?)
%  
% -------------------------------------------------------------------------
function [data, pseudoRollSmooth, swap_idx, check_idx] = ...
    estimatePseudoRoll(data, rotM_YP, largePertFlag, debugFlag, smoothType)
% -------------------------------------
%% inputs and params
if ~exist('rotM_YP','var') || isempty(rotM_YP)
    rotM_YP = [] ;
end
if ~exist('largePertFlag','var') || isempty(largePertFlag)
    largePertFlag = false ;
end
if ~exist('debugFlag','var') || isempty(debugFlag)
    debugFlag = false ;
end
if ~exist('smoothType','var') || isempty(smoothType)
    smoothType = 'sgolay' ; % 'sgolay' | 'lowpass' | 'wavelet'
end
% --------------
% params:
thetaB0_rad = 0 ; % rotate so that longitudinal body axis aligns with x axis
rollFlag = false ; % don't try to unroll the fly
frames = (data.params.startTrackingTime : data.params.endTrackingTime) ;
t = (1/data.params.fps)*frames ;
N_frames = length(frames) ;

diffThresh = 0.5*pi ; % 0.5*pi, 0.75*pi
minChunkLen = 15 ; % minimum length of data vector between two large derivative spikes that we'll consider
checkWindow = 6 ; % minimum length of data vector between two large derivative spikes that we'll consider
chunkDiffThresh = 0.70*pi ; 
zThresh = 2.5 ; % threshold for z-scored difference between smoothed and raw signals

% ------------------
% smoothing
sgolay_order = 2 ; % savitzky golay smoothing parameters
sgolay_len = 101 ; % 101

hampel_k = 11 ;
hampel_sigma = 3 ;

roll_filt_lvl = 100 ; % Hz
% --------------------------------------------------------------
%% center bodyCM at [0,0,0] and undo yaw and pitch rotations
if isempty(rotM_YP)
    % if we're not supplied with rotation matrices, jsut transform the data
    % structure
    data_bodyFrame = labToBodyFrame(data, largePertFlag,thetaB0_rad, ...
        rollFlag) ;
    leftWingCM = data_bodyFrame.leftWingCM ;
    rightWingCM = data_bodyFrame.rightWingCM ; 
    
    % grab the unpitch/unyaw matrices for later:
    rotM_YP = data_bodyFrame.bodyFrameRotMats ;
else
   % if we're supplied with rotation matrices, apply them (the benefit of
   %  this is that it saves us from having to re-calculate)
   leftWingCM = nan(N_frames, 3) ; 
   rightWingCM = nan(N_frames, 3) ; 
   for i = 1:N_frames
       % subtract off body CM and rotate
       R = squeeze(rotM_YP(:,:,i)) ; 
       leftWingCM(i,:) = R*(data.leftWingCM(i,:) - data.bodyCM(i,:))' ; 
       rightWingCM(i,:) = R*(data.rightWingCM(i,:) - data.bodyCM(i,:))' ; 
   end
end
    
% --------------------------------------------------------------
%% get vector from right to left wing
right2LeftVec = leftWingCM - rightWingCM ;
right2LeftVec = right2LeftVec./myNorm(right2LeftVec) ; % normalize

% project onto yz plane (perp. to body axis)
right2LeftVecProj = [zeros(size(right2LeftVec,1),1), right2LeftVec(:,2:3)] ;
right2LeftVecProj = right2LeftVecProj ./ myNorm(right2LeftVecProj) ;


% ---------------------------------------------------------------
%% estimate pseudo-roll from angle made by R2L vec
pseudoRollEst = atan2(right2LeftVecProj(:,3), right2LeftVecProj(:,2)) ;

% interpolate through nan values
nan_idx = isnan(pseudoRollEst) ;
pseudoRollEst = interp1(t(~nan_idx), pseudoRollEst(~nan_idx), t,'nearest') ;

% ------------------------------------
% hampel filter to detect outliers
[pseudoRollHampel, hampel_idx] = hampel(pseudoRollEst, hampel_k, ...
    hampel_sigma) ;
%pseudoRollHampel = pseudoRollEst ;

if (0)
    figure ;
    hold on
    plot(t, pseudoRollEst, 'ko')
    plot(t(hampel_idx), pseudoRollEst(hampel_idx),'rx')
    plot(t, pseudoRollHampel, 'b-')
    axis tight
end

% -----------------------------------------------------------------------
%% locate large spikes in derivative (correspond to potential L<->R swaps)
pseudoRollDiff = [0, diff(pseudoRollHampel)];
large_diff_idx = (abs(pseudoRollDiff) > diffThresh) ;

% get chunks of data in between large derivative spikes
large_diff_idx(1) = true ;
large_diff_idx(end) = true ;
chunk_ind_list = idx_by_thresh(~large_diff_idx) ;

% get length of each data chunk and exclude short ones (likely to be
% erroneous)
chunk_lengths = cellfun(@(y) length(y), chunk_ind_list) ;
too_short_idx = (chunk_lengths < minChunkLen) ;
short_chunk_ind = chunk_ind_list(too_short_idx) ; % save these for later so we can interpolate through
chunk_ind_list = chunk_ind_list(~too_short_idx) ; % remove short chunks for future calculations

data_chunk_list = cell(size(chunk_ind_list)) ;
for j = 1:length(chunk_ind_list)
    data_chunk_list{j} = pseudoRollHampel(chunk_ind_list{j}) ;
end

% ---------------------------------------------------------------------
%% check the initial chunk for pi shifts
% currently going by the logic that, in the body frame, the z coordinates
% of the wings should on average be above the body CM if the fly is right
% side up
ind_init = chunk_ind_list{1} ; 
wingCM_z = nanmean([rightWingCM(ind_init,3), leftWingCM(ind_init,3)],2) ; % average over left and right wing
wingCM_z_mean = nanmean(wingCM_z) ; % average over initial chunk
chunk1_roll_mean = nanmean(pseudoRollHampel(ind_init)) ;  

% if the mean wing CM z coordinate is less than zero, the fly is likely
% flipped over (rho ~ pi). If our estimate of roll thinks that the roll is
% zero in this interval, that indicates a mismatch that we need to correct,
% since later pi shifts rely on the previous chunks.
if (wingCM_z_mean < -4) && (abs(chunk1_roll_mean) < diffThresh)
    % if fly is flipped, we need to figure out whether to call that pi or
    % -pi. The best way to do this is probably using information from after
    % we've performed the other shifts by pi (then we'll just do a global
    % shift). So just flag it for now
    badFirstChunkFlag = true ;  
elseif (wingCM_z_mean > -4) && (abs(chunk1_roll_mean) > diffThresh)
    % opposite of above situation
    badFirstChunkFlag = true ;  
else
    badFirstChunkFlag = false ;  
end
% ------------------------------------------------------------------------
%% piece together new roll estimate from chunks
% new copy of data
pseudoRollInterp = pseudoRollHampel ; 

% ----------------------------------------
% shift chunks by +/- pi if the difference between the previous chunk's end
% and the current chunk's beginning is sufficiently large
%   NB: this assumes the first chunk is okay--need to check that
for k = 2:length(chunk_ind_list)
    ind_curr = chunk_ind_list{k} ; % current chunk indices
    ind_prev = chunk_ind_list{k-1} ; % previous chunk indices
    
    % read in roll estimates over these two chunks 
    data_prev = pseudoRollInterp(ind_prev) ;
    data_curr = pseudoRollInterp(ind_curr) ;
    
    % filter outliers that may make it difficult to sew together chunks
    data_prev = hampel(data_prev, minChunkLen) ;
    data_curr = hampel(data_curr, minChunkLen) ;
    
    % get (robust) estimates for the values of roll estimates at the
    % junction point between the chunks
    prev_right_end = nanmedian(data_prev((end-checkWindow+1):end)) ;
    curr_left_end = nanmedian(data_curr(1:checkWindow)) ;
    
    % if there is a jump between the previous and current jumps that is
    % order pi or greater, add or subtract multiples of pi until they fit
    % together
    while (curr_left_end - prev_right_end) > chunkDiffThresh
        data_curr = data_curr - pi ;
        curr_left_end = nanmedian(data_curr(1:checkWindow)) ;
    end
    while (curr_left_end - prev_right_end) < -1*chunkDiffThresh
        data_curr = data_curr + pi ;
        curr_left_end = nanmedian(data_curr(1:checkWindow)) ;
    end
    
    pseudoRollInterp(ind_curr) = data_curr ;
    
end

% assign nan values to both 1) short chunks and 2) points of large
% derivative
for kk = 1:length(short_chunk_ind)
    pseudoRollInterp(short_chunk_ind{kk}) = nan ;
end
pseudoRollInterp(large_diff_idx) = nan ;

% interpolate through remaining nan values
pseudoRollInterp = interp1(t(~isnan(pseudoRollInterp)), ...
    pseudoRollInterp(~isnan(pseudoRollInterp)), t, 'nearest','extrap') ;

% because idx_by_thresh sometimes misses some points, perform another
% hampel filter
pseudoRollInterp = hampel(pseudoRollInterp, hampel_k, hampel_sigma) ;

% finally, perform global shift if the first chunk seemed off
if badFirstChunkFlag  
    % should probably improve this check. for flies doing lots of rolling,
    % probably not always true that the median roll should be in the [-pi,
    % pi] range. but this should work for now
    pseudoRollMedian = nanmedian(pseudoRollInterp) ;
    if (pseudoRollMedian > diffThresh) 
        pseudoRollInterp = pseudoRollInterp - pi ; 
    elseif (pseudoRollMedian < -1*diffThresh) 
        pseudoRollInterp = pseudoRollInterp + pi ; 
    else
        fprintf('Cannot determine appropriate global shift \n')
        keyboard
    end
end
% ---------------------------------------------------------------------
%% plot results of pi shifts?
if debugFlag
    figure ;
    hold on
    plot(t, pseudoRollDiff)
    plot(t(large_diff_idx), pseudoRollDiff(large_diff_idx),'rx')
    axis tight
    
    figure ;
    hold on
    large_diff_ind = find(large_diff_idx) ;
    plot(t, pseudoRollHampel)
    axis tight
    ylim = get(gca,'ylim') ;
    for i = 1:length(large_diff_ind)
        plot(t(large_diff_ind(i)).*[1, 1], ylim, 'r--')
    end
    
    figure ;
    hold on
    for k = 1:length(chunk_ind_list)
        plot(t(chunk_ind_list{k}), pseudoRollHampel(chunk_ind_list{k}),'o')
    end
    %plot([t(1), t(end)], global_median*[1,1],'k--')
    plot(t, pseudoRollInterp,'r-')
    axis tight
end

% ---------------------------------------------------------
%% get smoothed signal
switch smoothType
    case 'wavelet'
        %wavelet denoising:
        pseudoRollSmooth = wdenoise(pseudoRollInterp) ;
        
    case 'lowpass'
        % lowpass filter:
        pseudoRollSmooth = filterEulerAngle(pseudoRollInterp, ...
            roll_filt_lvl ) ;
        
    case 'sgolay'
        % savitzky-golay:
        pseudoRollSmooth = sgolayfilt(pseudoRollInterp, sgolay_order, ...
            sgolay_len) ;
        
    otherwise
        fprintf('Invalid smooth type: %s \n', smoothType)
        keyboard
end

% ---------------------------------------------------------------
%% find points that significantly deviate from smoothed curve
% these will be the points where we either 1) swap wings, 2) add to ignore
% frames, or 3) just leave be
smoothDiff = pseudoRollSmooth - pseudoRollEst  ;
smoothDiffByPi = smoothDiff./pi ; % divide differences by pi
smoothDiffByPiInt = round(smoothDiffByPi) ; % round to nearest integer

% find data indices where number of pi rotations is odd. These need
% swapping
modPiRotations = mod(abs(smoothDiffByPiInt), 2) ;
swap_idx = (modPiRotations == 1) ; 

% next check for smaller differences between raw and smoothed roll--these
% may correspond to poorly tracked wings
pseudoRollSwap = pseudoRollEst + pi.*smoothDiffByPiInt ;
smoothDiffWithSwap = pseudoRollSmooth - pseudoRollSwap ; 

% take z score (or modified z score w median and MAD vs mean and STD)
% diffMedian = nanmedian(smoothDiffWithSwap) ; 
% diffMAD = nanmedian(abs(smoothDiffWithSwap - diffMedian)) ; 
% smoothDiff_z = (smoothDiffWithSwap - diffMedian)./ (diffMAD*1.4826) ; % the 1.4826 comes sigma ~ 1.4826*MAD (for norm dist)
smoothDiff_z = (smoothDiffWithSwap - nanmean(smoothDiffWithSwap))./ ...
    nanstd(smoothDiffWithSwap);

check_idx = (abs(smoothDiff_z) > zThresh) ;  

if (0)
    figure ;
    subplot(2,1,1)
    hold on
    plot(t, pseudoRollEst,'k.')
    plot(t, pseudoRollSwap)
    plot(t, pseudoRollSmooth)
    plot(t(check_idx), pseudoRollSwap(check_idx), 'rx')
    axis tight
    
    subplot(2,1,2)
    histogram(smoothDiff_z, 50)
end

% ----------------------------------------------------------------
%% swap wings that we flagged and add any new frames to ignore
swapFrames = find(swap_idx) ; 
data = swapWingLeftRight(data,swapFrames) ; 

ignoreFrames = find(check_idx) ; 
if isfield(data, 'ignoreFrames') && ~isempty(ignoreFrames)
    data.ignoreFrames = sort(unique([data.ignoreFrames, ignoreFrames])) ; 
elseif ~isfield(data, 'ignoreFrames') && ~isempty(ignoreFrames)
    data.ignoreFrames = ignoreFrames ; 
end

% ----------------------------------------------------------------
%% get new rollVectors and rhoTimes
rhoTimes = 1:N_frames ; 
rollVectors = repmat([0, 1, 0], N_frames, 1) ; 
% to get new roll vectors, we first define them in the body frame (by 
% rotating yhat in the yz plane), then in the lab frame, by re-applying
% body yaw and pitch rotations

for j = 1:N_frames 
    % value of roll angle from smoothed estimate
   rho = pseudoRollSmooth(j) ; 
   % matrices to first roll the vector, then rotate to lab frame
   rollRot = eulerRotationMatrix(0,0,rho)' ; 
   YPRot = squeeze(rotM_YP(:,:,j))' ; 
   
   % apply rotations
   rollVectors(j,:) = YPRot*rollRot*rollVectors(j,:)' ; 
end

% add updated values to struct
data.rhoTimes = rhoTimes ; 
data.rollVectors = rollVectors ; 

% ---------------------------------------------------------------
%% plot estimated roll?
if debugFlag
    % load previously estimated roll, if it exists (for visualization)
    defineConstantsScript
    try
        rhoRad = (pi/180)*data.anglesLabFrame(:,RHO) ;
    catch
        rhoRad = nan(size(t)) ;
    end
    
    % make plot of time series
    figure ;
    hold on
    plot(t, pseudoRollInterp, 'k.')
    plot(t, rhoRad,'b-')
    plot(t, pseudoRollSmooth,'r-')
    
    axis tight
    
end


end