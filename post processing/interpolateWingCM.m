%--------------------------------------------------------------------------
% function to interpolate through frames with bad estimates of wing cm.
% this is to be used for "cleanUpWingVoxels.m" where we try to correct for
% blobby or generally bad wing reconstructions
%
% INPUTS:
%       -data_in = standard analysis data structure
%       -wing_side = string that's either 'right' or 'left'
%       -debugFlag = boolean on whether or not to produce plots to check
%       interpolation
%--------------------------------------------------------------------------
function [wingCM_interp, good_idx, wingTip_interp] = ...
    interpolateWingCM(data, wing_side, bad_frames_init, debugFlag)
%--------------------------------------------------------------------------
% params and inputs
if ~exist('bad_frames_init','var') || isempty(bad_frames_init)
    bad_frames_init = false(data.Nimages,1) ; 
end
if ~exist('debugFlag','var') || isempty(debugFlag)
    debugFlag = false ; 
end

smoothType = 'rloess' ; 
smoothWindow = 7 ; %11 %9
diffThreshMin = 2 ; % vox/frame. just so if the median + MAD is too small, 
                    % we're not focusing on okay frames
minIndChunkSize = 3 ; % when we're looking at remaining good indices, how large of chunks to keep
hampelWindow = 71 ; 
hampelSigThresh = 2.5 ; 
interpDiffThresh = 5 ; 

if strcmpi(wing_side,'R')
    wing_side = 'right' ; 
elseif strcmpi(wing_side,'L')
    wing_side = 'left' ; 
end
% --------------------------------------------------------------------------
%% use jumps in the wing center of mass as an indicator of bad frames. 
% so if the voxel reconstruction for a single wing gives two blobs, this 
% should ideally show up as a blip in the movement

% initialize bad_ind (frames that we want to ignore for interpolation)
%bad_ind = false(data.Nimages,1) ; 

% load in data
wingCM = data.([wing_side 'WingCM']) ; 
wingTip = data.([wing_side 'WingTips']) ; 
bodyCM = data.bodyCM ; 
N_frames = data.Nimages ; 
frames = (1:N_frames)' ; 

t_frames = data.params.startTrackingTime : data.params.endTrackingTime ; 
zero_ind = find(t_frames > 0, 1, 'first') ; 

% subtract off body center of mass from wing cm and tip
wingCM = wingCM - bodyCM ;
wingTip = wingTip - bodyCM ; 

% exclude frames that we know to be bad
wingCM(bad_frames_init,:) = nan ; 

% -----------------------------------------------------------
%% first exclude points with high velocity
for i = 1:3
    % loop through x,y,z
    wingCMDiff = diff(wingCM(:,i)) ; 
    meanVel = nanmean(wingCMDiff(1:zero_ind)) ; 
    stdVel = nanstd(wingCMDiff(1:zero_ind)) ; 
    velOutlierIdx = [false; (abs(wingCMDiff) > (meanVel + 3*stdVel)) ]; 
    wingCM(velOutlierIdx,i) = NaN ; 
    
    % now interpolate
    %wingCM(:,i) = interp_nans(wingCM(:,i)) ; 
end
% find CM vec norm and interpolate nan values
normWingCM = myNorm(wingCM) ; 
% sometimes we remove too many points...
if (length(normWingCM) - sum(isnan(normWingCM))) < 10
    wingCM = data.([wing_side 'WingCM']) - bodyCM ;  
    normWingCM = myNorm(wingCM) ;
end

% --------------------------------------------------------------------
%% hampel filter for final outlier removal
[~, hampel_idx] = hampel(normWingCM, hampelWindow, hampelSigThresh) ; 
if (0)
    figure ; 
    plot(frames, normWingCM, 'o-')
    hold on 
    plot(frames(hampel_idx), normWingCM(hampel_idx),'rx')
end
normWingCM(hampel_idx) = nan ; 

%bad_ind = bad_ind | isnan(normWingCM) ; 
normWingCM = interp_nans(normWingCM) ; 

% smooth the CM vec norms using the 'rloess' option, which should be
% robust to outliers
normWingCM_smooth = smooth(normWingCM, smoothWindow, smoothType) ;

% compare smoothed and raw results--bad frames should have a large
% difference between the two
smoothDiff = abs(normWingCM - normWingCM_smooth) ; 

% ... to determine how large the difference needs to be for the frame to be
% considered bad, try to locate outliers using median + 3*MAD
smoothDiff_median = nanmedian(smoothDiff) ; 
smoothDiff_MAD = nanmedian(abs(smoothDiff - smoothDiff_median)) ;
smoothDiff_thresh = smoothDiff_median + 5*smoothDiff_MAD ; 

bad_ind = (smoothDiff > max([smoothDiff_thresh, diffThreshMin])) | ...
    bad_frames_init ;

% -------------------------------------------------------------------------
%% fit interpolant to good frames, to get predicted wing motion
% first find indices for frames that we'll use to fit interpolant
good_ind_list = idx_by_thresh(~bad_ind) ; 
good_ind_lengths = cellfun(@(y) length(y), good_ind_list) ; 
good_ind_list = good_ind_list(good_ind_lengths > minIndChunkSize) ;
good_ind = false(N_frames,1) ; 
for i = 1:length(good_ind_list)
   good_ind(good_ind_list{i}) = true ;  
end
good_ind = good_ind & ~any(isnan(wingCM),2) ; 

%--------------------------------------------------------
% now fit interpolants
wingCM_interp = nan(N_frames, 3) ; 
wingTip_interp = nan(N_frames, 3) ; 
for j = 1:3 
    % first fit interpolants
    c_cm = fit(frames(good_ind), wingCM(good_ind,j),'pchipinterp') ;
    wingCM_interp(:,j) = c_cm(frames) ; 
    nan_ind_tip = isnan(wingTip(:,j)) ; 
    c_tip = fit(frames(good_ind & ~nan_ind_tip), ...
        wingTip(good_ind & ~nan_ind_tip,j),'pchipinterp') ; 
    wingTip_interp(:,j) = c_tip(frames) ; 
    
    % then smooth
    wingCM_interp(:,j) = smooth(wingCM_interp(:,j), smoothWindow, ...
        smoothType) ;
    wingTip_interp(:,j) = smooth(wingTip_interp(:,j), smoothWindow, ...
        smoothType) ;
end
 
% -------------------------------------------------------------------------
%% plot some results?
if debugFlag
    %--------------------------------
    % time series of raw and smoothed
    figure ; 
    hold on
    plot(frames, normWingCM, 'o-')
    plot(frames, normWingCM_smooth, 'r-')
     plot(frames(bad_ind), normWingCM(bad_ind), 'ko',...
        'MarkerFaceColor','k')
    title([wing_side ' Wing CM Norm'])
    
    
    %--------------------------------
    % histogram of differences
%     figure ; 
%     nbins = 100 ;
%     histogram(smoothDiff, nbins)
%     %set(gca,'yscale','log')
%     title([wing_side ' Wing Smooth Diff'])
    
    %--------------------------------
    %interpolation check
    figure ; 
    for k = 1:3
        subplot(1,3,k)
        hold on
        plot(frames, wingCM(:,k), 'o-')
        plot(frames, wingCM_interp(:,k), 'r-')
        title([wing_side 'Wing Interpolation'])
    end
end

% --------------------------------------------------------
%% find good indices (ones we don't need to check)
% good indices are when the difference between raw and interpolated values
% is "small enough"

% add back in body center of mass
wingCM_interp = wingCM_interp + bodyCM ; 
wingTip_interp = wingTip_interp + bodyCM ; 

CM_interp_diff = myNorm(wingCM_interp - data.([wing_side 'WingCM'])) ;
good_idx = (CM_interp_diff < interpDiffThresh ) ; 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper functions
function y_out = interp_nans(y_in)
    x = (1:length(y_in))' ; 
    nan_idx = isnan(y_in) ; 
    y_out = interp1(x(~nan_idx), y_in(~nan_idx), x, 'linear', 'extrap') ; 
end