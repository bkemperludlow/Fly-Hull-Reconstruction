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
function [wingCM_interp, good_ind, wingTip_interp] = ...
    interpolateWingCM(data, wing_side, debugFlag)
%--------------------------------------------------------------------------
% params and inputs
if ~exist('debugFlag','var')
    debugFlag = false ; 
end

smoothType = 'rloess' ; 
smoothWindow = 7 ; %11 %9
diffThreshMin = 2 ; % vox/frame. just so if the median + MAD is too small, 
                    % we're not focusing on okay frames
minIndChunkSize = 5 ; % when we're looking at remaining good indices, how large of chunks to keep
% --------------------------------------------------------------------------
%% use jumps in the wing center of mass as an indicator of bad frames. 
% so if the voxel reconstruction for a single wing gives two blobs, this 
% should ideally show up as a blip in the movement

% load in data
wingCM = data.([wing_side 'WingCM']) ; 
wingTip = data.([wing_side 'WingTips']) ; 
N_frames = data.Nimages ; 
frames = (1:N_frames)' ; 

% find CM vec norm and interpolate nan values
normWingCM = myNorm(wingCM) ; 
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

bad_ind = (smoothDiff > max([smoothDiff_thresh, diffThreshMin])) ;

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
good_ind = good_ind & ~isnan(myNorm(wingCM)) ; 

%--------------------------------------------------------
% now fit interpolants
wingCM_interp = nan(N_frames, 3) ; 
wingTip_interp = nan(N_frames, 3) ; 
for j = 1:3 
    c_cm = fit(frames(good_ind), wingCM(good_ind,j),'smoothingspline') ; 
    wingCM_interp(:,j) = c_cm(frames) ; 
    c_tip = fit(frames(good_ind), wingTip(good_ind,j),'smoothingspline') ; 
    wingTip_interp(:,j) = c_tip(frames) ; 
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
    figure ; 
    nbins = 100 ;
    histogram(smoothDiff, nbins)
    %set(gca,'yscale','log')
    title([wing_side ' Wing Smooth Diff'])
    
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

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper functions
function y_out = interp_nans(y_in)
    x = (1:length(y_in))' ; 
    nan_idx = isnan(y_in) ; 
    y_out = interp1(x(~nan_idx), y_in(~nan_idx), x) ; 
end