%--------------------------------------------------------------------------
% function to correct errors in wing voxel reconstruction due to occlusion,
% etc. I'm beginning this to deal with cases in which the wings are large
% blobs because e.g. the side view of one wing is being counted as the side
% view of another.
%
% Note: currently, this will only work if most of the frames are okay and
% it's just e.g. the back of the wing stroke that has problems
%--------------------------------------------------------------------------
function data_out = cleanUpWingVoxels(data_in, all_fly_bw, ...
    body_only_bw, dlt, order, debugFlag)
%--------------------------------------------------------------------------
%{
 data_in = importdata(['D:\Fly Data\VNC Motor' ...
   'Lines\58_23032019\Analysis\Roll' ...
   'Left\Expr_58_mov_007\Expr_58_mov_007_test.mat']) ; 
%}
%% params and inputs
if ~exist('debugFlag','var')
    debugFlag = false ;
    figPosAfter = [2841, 273, 725, 610] ;
    figPosBefore = [ 2000, 270, 757, 613] ;
elseif exist('debugFlag','var') && debugFlag 
    figPosAfter = [2841, 273, 725, 610] ;
    figPosBefore = [ 2000, 270, 757, 613] ;
end

clustDebugFlag1 = false ;
clustDebugFlag2 = false ;

voxZscoreThresh = 1.5 ;
distThreshCM = 5 ; % vox
solidityThresh = 0.70 ; 
swapThresh = 10 ; % vox distance units
voxHogThresh = 0.80 ; % fraction of voxels taken by a single wing to merit clustering
voxProjThresh = 0.75 ; % minimum fraction of pixels that need to be hit for good wing assignment

% initialize output as copy of input
data_out = data_in ;
data_out.ignoreFrames = [] ; 
% i messed up with assignment of alternative chords--make sure to
% initialize them here
data_out.rightChordAltHats = data_out.chord1AltHats ; 
data_out.leftChordAltHats = data_out.chord2AltHats ; 
% -------------------------------------------------------------------------
%% load in some relevant data
N_frames = data_in.Nimages ;
rightWingInd = data_in.rightWingInd ;
leftWingInd = data_in.leftWingInd ;
bodyCM = data_in.bodyCM ;
wingLength = data_in.wingLength ;

params = data_in.params ; 
voxelSize = params.voxelSize ; 
detectorLengthPix = params.detectorLengthPix ;
CAMERAS = params.CAMERAS ;
NCAMS = params.NCAMS ;

% get row indices
df = diff(data_in.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data_in.res,1)] ;
clear df ;
% -------------------------------------------------------------------------
%% get voxel counts to try to diagnose issues with frames
N_vox_R = nan(N_frames, 1) ;
N_vox_L = nan(N_frames, 1) ;
for j = 1:N_frames
    row1 = frameStartInd(j) ;
    row2 = frameEndInd(j) ;
    
    IDX = data_in.RESIDX(row1:row2,:) ; % Classifying each voxel as L/R wing or body
    wingRows_R = (IDX(:,rightWingInd)==1) ;
    wingRows_L = (IDX(:,leftWingInd)==1) ;
    N_vox_R(j) = sum(wingRows_R) ;
    N_vox_L(j) = sum(wingRows_L) ;
end

% z-score the voxel numbers
mean_vox_R = nanmean(N_vox_R) ;
std_vox_R = nanstd(N_vox_R) ;
mean_vox_L = nanmean(N_vox_L) ;
std_vox_L = nanstd(N_vox_L) ;

N_vox_R_zScore = (N_vox_R - mean_vox_R)./std_vox_R ;
N_vox_L_zScore = (N_vox_L - mean_vox_L)./std_vox_L ;

% get fraction of total wing voxels per frame that are assigned to left vs
% right
N_vox_R_frac = N_vox_R./(N_vox_R + N_vox_L) ; 
N_vox_L_frac = 1 - N_vox_R_frac ; 

% tests to determine if we should expect merged wings
rightVoxHogFlag = (N_vox_R_frac > voxHogThresh) ;
leftVoxHogFlag = (N_vox_L_frac > voxHogThresh) ;
voxHogFlag = rightVoxHogFlag | leftVoxHogFlag ; % is one wing hogging all the voxels?
enoughVoxFlag = (N_vox_L_zScore > 0) | (N_vox_R_zScore > 0) ; % are there enough wing voxels to merit a split?
guessMergeFlag = voxHogFlag & enoughVoxFlag ;

% -------------------------------------------------------------------------
%% use wing center of mass jumps to determine putative bad recon. frames
[rightWingCM_interp, good_idx_R, rightWingTips_interp] = ...
    interpolateWingCM(data_in,'right', guessMergeFlag) ;
[leftWingCM_interp, good_idx_L, leftWingTips_interp] = ...
    interpolateWingCM(data_in,'left', guessMergeFlag) ;

rightWingTips_interp = filter_wing_vecs(rightWingTips_interp,3,2,2000) ;
leftWingTips_interp = filter_wing_vecs(leftWingTips_interp,3,2,2000) ;

bad_idx = ~good_idx_R | ~good_idx_L ;
bad_frames = find(bad_idx) ;
% can't correct first frame because need info from previous frames
if (sum(bad_frames == 1) > 0)
    bad_frames = bad_frames(2:end) ; 
    bad_idx(1) = false ; 
    data_out.ignoreFrames = [data_out.ignoreFrames, 1] ; 
end
N_bad_frames = length(bad_frames) ;

% -------------------------------------------------------------------------
%% loop through bad frames and see if we can sort out the issue
if debugFlag
    h_debug_in = figure('PaperPositionMode','auto','Position',figPosBefore) ;
    h_debug_out = figure('PaperPositionMode','auto','Position',figPosAfter) ;
end

clusterRightFlag = nan(N_bad_frames,1) ;
clusterLeftFlag = nan(N_bad_frames,1) ;
clusterErrorFlag = nan(N_bad_frames, 1) ;

for k = 1:N_bad_frames 
    ind = bad_frames(k) ;
    fprintf('Frame %d : beginning analysis \n', ind )
    row_start = frameStartInd(ind) ; % rows for voxel idx array
    row_end = frameEndInd(ind) ;
    %----------------------------------------------------------------------
    % generate plots to check situation?
    if debugFlag
        h_debug_in = plotFlyVoxels(data_in, ind, h_debug_in) ;
        title(['Frame ' num2str(ind) ' BEFORE'])
    end
    
    %----------------------------------------------------------------------
    %% get segmented wing images from all views
    % first need wing voxels to assign right vs left
    [wingVoxR, wingVoxL] = readWingVox(data_in, row_start, row_end) ;
    
    % segment wings
    [imWingRMat, imWingLMat, segWingsFlag] = ...
        getSegmentedWingImgs(ind, all_fly_bw,  body_only_bw, wingVoxR, ...
        wingVoxL,  dlt, params, order) ;
    
    % ---------------------------------------------------------------------
    %% check if the same voxels have been assigned to both wings
    % sometimes a weird thing happens where both the right and left wings
    % are assigned the same voxels, so they look like they overlap
    numUniqueVox = size(unique([wingVoxR ; wingVoxL],'rows'),1) ; 
    voxOverlapCheck = (numUniqueVox == size(wingVoxR,1)) & ...
        (numUniqueVox == size(wingVoxL,1)) & (numUniqueVox > 0);
    if voxOverlapCheck
       % if this happens, the same voxels are likely assigned to both
       % wings. need to figure out which wing the voxels should belong to
       fprintf('Frame %d : same voxels assigned to both wings \n', ind )
       wingVox = unique([wingVoxR ; wingVoxL],'rows') ;
       if (sum(segWingsFlag) > 0)
           projR = nan(1,NCAMS) ;
           projL = nan(1,NCAMS) ;
           for camNum = CAMERAS(segWingsFlag)
               % grab wing images
               imWingR = squeeze(imWingRMat(camNum,:,:)) ;
               imWingL = squeeze(imWingLMat(camNum,:,:)) ;
               
               % project vox onto pix
               projR(camNum) = checkVoxAgainstImg(wingVox, imWingR,...
                   dlt, camNum, detectorLengthPix, order, voxelSize) ;
               projL(camNum) = checkVoxAgainstImg(wingVox, imWingL,...
                   dlt, camNum, detectorLengthPix, order, voxelSize) ;
           end
           
           % compare projections across cameras
           meanProjR = nanmean(projR) ; 
           meanProjL = nanmean(projL) ; 
           
           % figure out which voxels to assign 
           if (meanProjR >= meanProjL)
               delete_wing = 'left' ; 
               N_vox_L_frac(ind) = 0 ; 
               N_vox_R_frac(ind) = 1 ; 
               N_vox_L_zScore(ind) = -Inf ; 
           elseif (meanProjR < meanProjL)
               delete_wing = 'right' ; 
               N_vox_R_frac(ind) = 0 ; 
               N_vox_L_frac(ind) = 1 ; 
               N_vox_R_zScore(ind) = -Inf ; 
           end
       else
           % if we don't have good segmented wing image, use centroid
           voxCentroid = nanmean(wingVox,1) ; 
           CM_distR = myNorm(rightWingCM_interp(ind,:) - voxCentroid) ; 
           CM_distL = myNorm(leftWingCM_interp(ind,:) - voxCentroid) ; 
           
           % figure out which voxels to assign 
           if (CM_distR <= CM_distL)
               delete_wing = 'left' ; 
               N_vox_L_frac(ind) = 0 ; 
               N_vox_R_frac(ind) = 1 ; 
               N_vox_L_zScore(ind) = -Inf ; 
           elseif (CM_distR > CM_distL)
               delete_wing = 'right' ; 
               N_vox_R_frac(ind) = 0 ; 
               N_vox_L_frac(ind) = 1 ; 
               N_vox_R_zScore(ind) = -Inf ; 
           end
       end
       
       % ------------------------------------------------------------------
       % now remove voxels/wing data from the unassigned wing. the assigned
       % wing should still be fine
       data_out.RESIDX(row_start:row_end, ...
           data_out.([delete_wing 'WingInd'])) = false ;
       data_out.([delete_wing 'SpanHats'])(ind,:) = nan(1,3) ; 
       data_out.([delete_wing 'ChordHats'])(ind,:) = nan(1,3) ;
       data_out.([delete_wing 'WingCM'])(ind,:) = nan(1,3) ;
       data_out.([delete_wing 'WingTips'])(ind,:) = nan(1,3) ;
       
       % since this means that all of the voxels belong to one wing, we
       % should try to unmerge the wing voxels, provided there's enough
       guessMergeFlag(ind) = (N_vox_L_zScore(ind) > 0) | ...
           (N_vox_R_zScore(ind) > 0) ;
    end
    %----------------------------------------------------------------------
    %% if wings are likely merged, cluster them
    % If all flags are tripped, we guess that 1) the wings
    % are merged and 2) there's a crappy little bundle of pixels for one
    % wing off in a corner. in this case, we'll cluster
    if guessMergeFlag(ind)
        fprintf('Frame %d : wings seemed merged, attempting to cluster... \n',...
            ind )
        if (N_vox_R_frac(ind) >= 0.5)
            wing_str = 'right' ;
        elseif (N_vox_R_frac(ind) < 0.5)
            wing_str = 'left' ;
        else
            disp('this should not happen')
            keyboard
        end
        % perform clustering
        [wingVox, wingRows, label_idx, centroids, badClusterFlag] = ...
            clusterWings(data_out, ind, wing_str, clustDebugFlag1 ) ; %clustDebugFlag1
       
        if badClusterFlag
            clusterErrorFlag(k) = true ;
            continue
        end
        
        % give the debug plot a title, if it comes up
        if clustDebugFlag1
            title(['Frame ' num2str(ind) ' merged wing clustering'])
        end
        
        % -------------------------------------------------------------
        % if we did successfully cluster, we now need to figure out
        % which is left vs right. we'll do this by projecting to binary 
        % images NB: old version used checking against the interpolated 
        % center of mass; code at bottom
        if (sum(segWingsFlag) > 0)
            N_clusts = size(centroids,1) ;
            projScoreR = nan(N_clusts, NCAMS) ;
            projScoreL = nan(N_clusts, NCAMS) ;
            for i = 1: N_clusts
                for camNum = CAMERAS(segWingsFlag)
                    % grab wing images
                    imWingR = squeeze(imWingRMat(camNum,:,:)) ;
                    imWingL = squeeze(imWingLMat(camNum,:,:)) ;
                    
                    % project vox onto pix
                    projScoreR(i,camNum) = ...
                        checkVoxAgainstImg(wingVox(label_idx == i,:), ...
                        imWingR, dlt, camNum, detectorLengthPix, order,...
                        voxelSize) ;
                    projScoreL(i,camNum) = ...
                        checkVoxAgainstImg(wingVox(label_idx == i,:), ...
                        imWingL, dlt, camNum, detectorLengthPix, order,...
                        voxelSize) ;
                    
                end
            end
            
            % get average projected fractions
%             meanProjR = nanmean(projScoreR, 2) ;
%             meanProjL = nanmean(projScoreL, 2) ;
            
            % since the wings are merged, we can't be super confident about
            % L/R assignments here. instead, just pick the view that gives 
            % the maximum projections
            [~, projCam] = nanmax( nanmax(projScoreR) + nanmax(projScoreL)) ;
            [rightProjMax, right_idx] = nanmax(projScoreR(:,projCam)) ;
            [leftProjMax, left_idx] = nanmax(projScoreL(:,projCam)) ;
%             rightProjMax = nanmax(projScoreR(right_idx,:)) ; 
%             leftProjMax = nanmax(projScoreL(left_idx,:)) ; 
            projThreshCheck = (rightProjMax >= voxProjThresh) && ...
                (leftProjMax >= voxProjThresh) ;
            
            unMergeFailFlag = ~(projThreshCheck && ~(right_idx == left_idx)) ;
        else
            % if we don't have any good segmented wing images
            dist_mat = pdist2(centroids, [rightWingCM_interp(ind,:) ; ...
                leftWingCM_interp(ind,:)]) ;
            diagSum = sum(diag(dist_mat)) ;
            offDiagSum = sum(diag(fliplr(dist_mat))) ;
            if (offDiagSum - diagSum) > swapThresh
                right_idx = 1 ;
                left_idx = 2 ;
                unMergeFailFlag = false ;
            elseif (offDiagSum - diagSum) < -1*swapThresh
                right_idx = 2 ;
                left_idx = 1 ;
                unMergeFailFlag = false ;
            else
                right_idx = 1 ;
                left_idx = 1 ;
                unMergeFailFlag = true ;
            end
        end
        
        % ----------------------------------------------------------
        % perform assignment:
        if unMergeFailFlag %|| ~new_dist_check
            % in case both blobs think they should be the same wing
            fprintf('Frame %d : clustering likely unnecessary or was inconclusive \n', ind )
        else
            % otherwise, assign the new wing data to our output struct
            wingVoxR = wingVox(label_idx == right_idx,:) ;
            rightCM = centroids(right_idx,:) ;
            wingRowsR = wingRows(:,right_idx) ;
            wingVoxL = wingVox(label_idx == left_idx,:) ;
            leftCM = centroids(left_idx,:) ;
            wingRowsL = wingRows(:,left_idx) ;
            
            % now calculate wing vectors and store them in data struct so
            % we can proceed with other checks
            rightRefVecs = [bodyCM(ind,:); bodyCM(ind-1,:); ...
                rightWingTips_interp(ind-1,:)] ;
            [spanHatR, chordHatR, chordAltHatR, ~, wingTipR] = ...
                estimate_wing_vecs(wingVoxR, rightRefVecs, wingLength, ...
                [], [], rightCM) ;
            
            leftRefVecs = [bodyCM(ind,:); bodyCM(ind-1,:); ...
                leftWingTips_interp(ind-1,:)] ;
            [spanHatL, chordHatL, chordAltHatL, ~, wingTipL] = ...
                estimate_wing_vecs(wingVoxL, leftRefVecs, wingLength, ...
                [], [], leftCM) ;
            
            % add to storage arrays
            data_out.rightWingCM(ind,:) = rightCM ;
            data_out.rightSpanHats(ind,:) = spanHatR ;
            data_out.rightChordHats(ind,:) = chordHatR ;
            data_out.rightChordAltHats(ind,:) = chordAltHatR ;
            data_out.rightWingTips(ind,:) = wingTipR ;
            N_vox_R_zScore(ind) = (size(wingVoxR,1) - mean_vox_R)/std_vox_R ;
            
            data_out.leftWingCM(ind,:) = leftCM ;
            data_out.leftSpanHats(ind,:) = spanHatL ;
            data_out.leftChordHats(ind,:) = chordHatL ;
            data_out.leftChordAltHats(ind,:) = chordAltHatL ;
            data_out.leftWingTips(ind,:) = wingTipL ;
            N_vox_L_zScore(ind) = (size(wingVoxL,1) - mean_vox_L)/std_vox_L ;
            
            %...including fucking voxels, blurg
            data_out.RESIDX(row_start:row_end,2:3) = [wingRowsR, wingRowsL] ;
            
            % make sure to store whether or not we clustered
            if strcmp(wing_str,'right')
                clusterRightFlag(k) = ~badClusterFlag ;
            elseif strcmp(wing_str,'left')
                clusterLeftFlag(k) = ~badClusterFlag ;
            end
        end
    else
        fprintf('Frame %d : no need to un-merge wings \n',ind)
    end
    
    %----------------------------------------------------------------------
    %% if a single wing is blobby try to cluster that too
    % ---------------------------------------------------------------
    %grab wing voxels again, so we include clustering changes above
    [wingVoxR, wingVoxL] = readWingVox(data_out, row_start, row_end) ;
    
    % ------------------------------------------------
    % perform some checks to see if we should cluster
    % ------------
    % right wing:
    if ~isempty(wingVoxR)
        checkRightCM = (myNorm(data_out.rightWingCM(ind,:) - ...
            rightWingCM_interp(ind,:)) > distThreshCM) ;
        justRightCheck = checkRightCM | (N_vox_R_zScore(ind) > voxZscoreThresh) ;

        [rightWingSolidity, ~, ~] = wingVoxSolidity(wingVoxR) ;
        solidityCheckR = (rightWingSolidity > solidityThresh) & ...
            (N_vox_R_zScore(ind) < voxZscoreThresh) ;
    else
        justRightCheck = false ; 
        solidityCheckR = false ;
    end
    clusterJustRightFlag = justRightCheck && ~solidityCheckR ; 
    
    % ------------
    % left wing:
    if ~isempty(wingVoxL)
        checkLeftCM = (myNorm(data_out.leftWingCM(ind,:) - ...
            leftWingCM_interp(ind,:)) > distThreshCM) ;
        justLeftCheck = checkLeftCM | (N_vox_L_zScore(ind) > voxZscoreThresh) ;
        
        [leftWingSolidity, ~, ~] = wingVoxSolidity(wingVoxL) ;
        solidityCheckL = (leftWingSolidity > solidityThresh) & ...
            (N_vox_L_zScore(ind) < voxZscoreThresh) ;
    else
        justLeftCheck = false ; 
        solidityCheckL = false ;
    end
    clusterJustLeftFlag = justLeftCheck && ~solidityCheckL ; 
    
    % --------------------------------------------------------------------
    % if at least one wing needs to be clustered, get images of segmented
    % wings to compare against candidate clusters
    % NB: we're redoing this because we may have clustered merged wings
    if clusterJustRightFlag || clusterJustLeftFlag
        [imWingRMat, imWingLMat, segWingsFlag] = ...
            getSegmentedWingImgs(ind, all_fly_bw,  body_only_bw, wingVoxR, ...
            wingVoxL,  dlt, params, order) ; 
    end
    oneWingClustFlag = false ;
    % ---------------------------
    % RIGHT wing
    if clusterJustRightFlag
        fprintf('Frame %d : clustering just right wing voxels \n', ind )
        wing_str = 'right' ;
        try
            [data_out, ~, ignoreFlagR] = clusterSingleWing(data_out, ...
                wing_str, ind, rightWingCM_interp, rightWingTips_interp, ...
                row_start, row_end, imWingRMat, segWingsFlag, dlt, order,...
                voxelSize, clustDebugFlag2) ;
        catch
            keyboard
        end
        oneWingClustFlag = true ;
        if ignoreFlagR
            data_out.ignoreFrames = unique([data_out.ignoreFrames, ind]) ;
        end
    end
    % ---------------------------
    % LEFT wing
    if clusterJustLeftFlag
        fprintf('Frame %d : clustering just left wing voxels \n', ind )
        wing_str = 'left' ;
        try
            [data_out, ~, ignoreFlagL] = clusterSingleWing(data_out, ...
                wing_str, ind, leftWingCM_interp, leftWingTips_interp, ...
                row_start, row_end, imWingLMat, segWingsFlag, dlt, order,...
                voxelSize, clustDebugFlag2) ;
        catch
            keyboard
        end
        oneWingClustFlag = true ;
        if ignoreFlagL
            data_out.ignoreFrames = unique([data_out.ignoreFrames, ind]) ;
        end
    end
    if ~oneWingClustFlag
        fprintf('Frame %d : no blobby wings to cluster \n', ind)
    end

    %----------------------------------------------------------------------
    %% generate plots to check our work?
    if debugFlag
        h_debug_out = plotFlyVoxels(data_out, ind, h_debug_out) ;
        title(['Frame ' num2str(ind) ' AFTER'])
        keyboard
    end
    fprintf('Completed %d /%d frames \n',k, N_bad_frames)
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
%% cluster one wing if it's too blobby
function [data_out, badClusterFlag, ignoreFlag] = ...
    clusterSingleWing(data_in, wing_side, ind, wingCM_interp, ...
    wingTip_interp, row_start, row_end, imWingMat, segWingsFlag, dlt, order,...
    voxelSize, clusterDebugFlag)
%--------------------------
% inputs and params
% data_out starts as a copy of data_in
data_out = data_in ;
ignoreFlag = false ; 

% load in a few params
wingInd = data_in.([wing_side 'WingInd']) ;
spanPrev = data_in.([wing_side 'SpanHats'])(ind,:) ; 
wingLength = data_in.wingLength ;
minNumVox = 500 ;
tipDistThresh = 2 ; 
pcaRatioThresh = 5 ; 
dotThresh = 0.9 ; 
%cmCompThresh = -2 ; 
%--------------------------------------------------------------------------
% use the "clusterWings" function to get the two best clusters from the vox
[wingVox, wingRows, label_idx, centroids, badClusterFlag] = ...
    clusterWings(data_out, ind, wing_side, clusterDebugFlag) ;

if clusterDebugFlag 
   title(['Frame ' num2str(ind) ' ' wing_side ' wing clustering']) 
end
if badClusterFlag
    fprintf('Frame %d : error performing clustering on %s wing \n', ind , wing_side)
    ignoreFlag = true ; 
    return
end

%--------------------------------------------------------------------------
% calculate vectors for putative wing clusters.
N_clusts = size(centroids,1) ;
spanHatMat = nan(N_clusts,3) ;
chordHatMat = nan(N_clusts,3) ;
chordAltHatMat = nan(N_clusts,3) ;
wingTipMat = nan(N_clusts,3) ;
pcaCoeffCell = cell(N_clusts,1) ; 
pcaLatentCell = cell(N_clusts,1) ; 
projScore = nan(N_clusts,3) ; 

refVecs = [data_in.bodyCM(ind,:); data_in.bodyCM(ind-1,:) ; ...
    wingTip_interp(ind-1,:)] ;

for i = 1:N_clusts
    [spanHat, chordHat, chordAltHat, ~, wingTip] = ...
        estimate_wing_vecs(wingVox(label_idx == i,:), refVecs, wingLength, ...
        [], [], centroids(i,:)) ;
    spanHatMat(i,:) = spanHat ;
    chordHatMat(i,:) = chordHat ;
    chordAltHatMat(i,:) = chordAltHat ;
    wingTipMat(i,:) = wingTip ;
    
    % also get pca info
    [coeff,~,latent] = pca(wingVox(label_idx == i,:)) ;
    pcaCoeffCell{i} = coeff ; 
    pcaLatentCell{i} = latent ; 
    
    % also check projection onto images
    for camNum = 1:3 
        if segWingsFlag(camNum)
            imWing = squeeze(imWingMat(camNum,:,:)) ; 
            detectorLengthPix = size(imWing,1) ; 
            projScore(i,camNum) = ...
                checkVoxAgainstImg(wingVox(label_idx == i,:), imWing, dlt, ...
                camNum, detectorLengthPix, order, voxelSize) ;
        end
    end
end

%------------------------------------------
% determine which cluster gives best wing
% NB: just taking the longest for now, should fix soon. maybe compare
% against XY view BW image?
tip2span_dists = myNorm(wingTipMat - centroids) ;
cm_dists = myNorm(repmat(wingCM_interp(ind,:),N_clusts,1) - centroids) ;
pca_ratio = cellfun(@(y) y(1)/y(end), pcaLatentCell) ; 
pca_span_dot = cellfun(@(y) abs(dot(spanPrev,y(:,1))), pcaCoeffCell) ;
mean_proj_frac = nanmean(projScore,2) ;
% pca_axis_cross_z = cellfun(@(y) y(3), pca_axis_cross) ;

[~, max_dist_idx] = max(tip2span_dists) ;
[~, closest_cm_idx] = min(cm_dists) ; 
[~, max_pca_idx] = max(pca_ratio) ;
[~, max_proj_idx] = max(mean_proj_frac) ; 
any_nan_idx = any(isnan(centroids),2) | any(isnan(spanHatMat),2) | ...
    any(isnan(chordHatMat),2) ; 
% if abs(tip2span_dists(1) - tip2span_dists(2)) > tipDistThresh
%     good_idx = max_dist_idx ; 
if sum(any_nan_idx) > 0
    tmp_idx = find(~any_nan_idx) ; 
    if (length(tmp_idx) == 1)
        good_idx = tmp_idx ; 
    else
        disp('Could not assign wing vectors, skipping...')
        ignoreFlag = true ; 
        return
    end
elseif (abs(mean_proj_frac(1) - mean_proj_frac(2)) > 0)
    good_idx = max_proj_idx ;
elseif (abs(pca_ratio(1) - pca_ratio(2)) > pcaRatioThresh) && ...
        (pca_span_dot(max_pca_idx) > dotThresh)
    good_idx = max_pca_idx ;
elseif abs(tip2span_dists(1) - tip2span_dists(2)) > tipDistThresh
     good_idx = max_dist_idx ; 
else
    good_idx = closest_cm_idx ; 
end
%-------------------------------------------------------------------
% check that the clustered bit even changes/improves the end result
% (for now, just test that there's a sufficient number of voxels and that
% we're not just splitting a wing down the middle, spanwise
% TO DO: should i just be comparing the xy projections? that's a prety
% reliable view, and splitting the wing spanwise would be flagged by that
interpSpan = wingTip_interp(ind,:) - wingCM_interp(ind,:) ; 
interpSpanHat = interpSpan ./ myNorm(interpSpan) ; 
centroidDispVec = centroids(1,:) - centroids(2,:) ; 
centroidDispHat = centroidDispVec ./ myNorm(centroidDispVec) ; 
dotCheck = dot(centroidDispHat, interpSpanHat) ; 
% cm_dist_comp = repmat(cm_dist_orig,N_clusts,1) - cm_dists ; 
if (sum(label_idx == good_idx) < minNumVox) || (abs(dotCheck) > dotThresh)
    disp('bad clustering--just take original')
else
    %----------------------------------
    % else enter data into structure
    data_out.([wing_side 'WingCM'])(ind,:) = centroids(good_idx,:) ;
    data_out.([wing_side 'SpanHats'])(ind,:) = spanHatMat(good_idx, :) ;
    data_out.([wing_side 'ChordHats'])(ind,:) = chordHatMat(good_idx,:) ;
    data_out.([wing_side 'ChordAltHats'])(ind,:) = chordAltHatMat(good_idx,:) ;
    data_out.([wing_side 'WingTips'])(ind,:) = wingTipMat(good_idx,:) ;
    data_out.RESIDX(row_start:row_end,wingInd) = wingRows(:,good_idx) ;
end

end

% ------------------------------------------------------------------------
%% check wing voxels against the binarized image of the wing
function projFrac = checkVoxAgainstImg(wingCoords, imWing, dlt, ...
    camNum, detectorLengthPix, order, voxelSize)
% -------------------------
% inputs
% if ~exist('detectorLengthPix','var') || isempty(detectorLengthPix)
%    detectorLengthPix = 512 ;  
% end
% if ~exist('order','var') || isempty(order)
%     order = [2, 1, 3] ; 
% end
% ---------------------------------------------------------------------
% make sure wing coordinates are in "real space" and not "voxel space"
if isinteger(wingCoords) || (max(wingCoords(:)) > 1e-2)
   wingCoords = voxelSize*double(wingCoords) ; 
end
% -----------------------------------------------------------
% project wing coordinates to image space with dlt
wingProjPix = dlt_inverse(dlt(:,order(camNum)), wingCoords) ; 
pixel_floor = floor([wingProjPix(:,1), ...
    (detectorLengthPix - wingProjPix(:,2) + 1)]) ; 
pixel_ceil = ceil([wingProjPix(:,1), ...
    (detectorLengthPix - wingProjPix(:,2) + 1)]) ;

% take out any that go outside image
good_idx_floor = (pixel_floor(:,1) > 0) & (pixel_floor(:,1) < detectorLengthPix) & ...
    (pixel_floor(:,2) > 0) & (pixel_floor(:,2) < detectorLengthPix) ;
good_idx_ceil = (pixel_ceil(:,1) > 0) & (pixel_ceil(:,1) < detectorLengthPix) & ...
    (pixel_ceil(:,2) > 0) & (pixel_ceil(:,2) < detectorLengthPix) ;

pixel_floor(~good_idx_floor,:) = 1 ; 
pixel_ceil(~good_idx_ceil,:) = 1 ; 

% get image indices for pixel values
pix_floor_idx = sub2ind(size(imWing), pixel_floor(:,2), pixel_floor(:,1)) ; 
pix_ceil_idx = sub2ind(size(imWing), pixel_ceil(:,2), pixel_ceil(:,1)) ;

% if we're checking this projection to trim voxels, do that. otherwise just
% check how good the projection is
imWingProj = false(size(imWing)) ; 
imWingProj(pix_floor_idx) = true ; 
imWingProj(pix_ceil_idx) = true ; 

% check image overlap
if (0)
   figure ; 
   imshowpair(imWing, imWingProj)
end
% comapare projection and binary image
projFrac = sum(sum(imWingProj & imWing))/ sum(imWing(:)) ; 

end


% ------------------------------------------------------------------------
%% compare different L/R assignments for a segmented wing images
% i.e. should "wing1" be right or left, ditto "wing2"
function dist_mat = compareWingCentroidVox(imWing1, imWing2, camNum, voxR,...
    voxL, dlt, order, voxelSize)
% initialize some storage
dist_mat = Inf*ones(2, 2) ; 
detectorLengthPix = size(imWing1,1) ; 

% concatenate arrays so we can loop through
imWingCell = {imWing1, imWing2} ; 
voxCell = {voxR, voxL} ; 
for ii = 1:2
    % ii corresponds to wing image number
    imWing = imWingCell{ii} ; 
    for jj = 1:2
        % jj corresponds to right vs left wing voxels
        wingVox = voxCell{jj} ; 
        
        % -----------------------------------------------------------------
        % project voxels onto camera image
        projFrac = checkVoxAgainstImg(wingVox, imWing, dlt, ...
            camNum, detectorLengthPix, order, voxelSize) ; 
        
        % calculate distance as one minus fraction of wing pixels hit by
        % projection
        dist_mat(ii,jj) = 1 - projFrac ; 
    end
end
end

% -------------------------------------------------------------------------
%% get segmented wing images
% putting this here to avoid code clutter 
function [imWingRMat, imWingLMat, segWingsFlag] = ...
    getSegmentedWingImgs(ind, all_fly_bw,  body_only_bw, wingVoxR, ...
        wingVoxL,  dlt, params, order)
% -----------------------------------
%% inputs and params
detectorLengthPix = params.detectorLengthPix ;
CAMERAS = params.CAMERAS ;
NCAMS = params.NCAMS ;
voxelSize = params.voxelSize ;
% -------------------------------------------------------------------------
% load in binary images of 1) full fly 2) just body
imMatFly   = false(NCAMS, detectorLengthPix, detectorLengthPix) ;
imMatBody  = false(NCAMS, detectorLengthPix, detectorLengthPix) ;
for c=CAMERAS
    % load images from sparse array
    imMatFly(c,:,:)  = getImage4D(all_fly_bw, c, ind);
    imMatBody(c,:,:) = getImage4D(body_only_bw, c, ind);
end

% segment wings in images
imWingRMat = false(NCAMS, detectorLengthPix, detectorLengthPix) ;
imWingLMat = false(NCAMS, detectorLengthPix, detectorLengthPix) ;
sameMasksFlag = false(NCAMS,1) ;
noWingsFlag = false(NCAMS,1) ;
wingExtremaFlag = false(NCAMS,2,2) ; % N_cams x N_objects x N_coordinates
% loop through cameras and try to segment
for c = CAMERAS
    imFly = squeeze(imMatFly(c,:,:)) ;
    imBody = squeeze(imMatBody(c,:,:)) ;
    [imWingRMat(c,:,:), imWingLMat(c,:,:), sameMasksFlag(c), noWingsFlag(c),...
        ~, ~, wingExtremaFlag(c,:,:)] = segmentWings_Sam(imFly, imBody ,...
            false) ;
end

% try to flag wings that are well-segmented
segWingsFlag = ~sameMasksFlag & ~noWingsFlag & ...
    all(any(wingExtremaFlag,3),2) ;

% initialize storage for projected fraction
% projFracR = nan(3,1) ; 
% projFracL = nan(3,1) ; 
% ------------------------------------------------------------
%% swap wings?
% I've just called wing1 the right wing and wing2 the left wing, but we
% need to actually check that
check_cams = CAMERAS(segWingsFlag) ;
for c = check_cams
    % similarity found using projection of voxels onto image
    imWingR = squeeze(imWingRMat(c,:,:)) ; 
    imWingL = squeeze(imWingLMat(c,:,:)) ; 
    dist_mat = compareWingCentroidVox(imWingR, imWingL, c, wingVoxR,...
        wingVoxL, dlt, order, voxelSize);
    % if the off diagonal elements have larger sum, this means there's
    % more overlap with "1" vs "2" switched in the side view
    if sum(diag(dist_mat)) > sum(diag(fliplr(dist_mat)))
        % switch images
        temp2 = squeeze(imWingRMat(c,:,:)) ;
        imWingRMat(c,:,:) = squeeze(imWingLMat(c,:,:)) ;
        imWingLMat(c,:,:) = temp2 ;
        %clear temp temp2
%         projFracR(c) = 1-dist_mat(2,1) ; 
%         projFracL(c) = 1-dist_mat(1,2) ; 
%     else
%         projFracR(c) = 1-dist_mat(1,1) ; 
%         projFracL(c) = 1-dist_mat(2,2) ; 
    end
end
end

% -------------------------------------------------------------------------
%% load wing voxels (just to declutter, should check if this slows code)
function [wingVoxR, wingVoxL] = readWingVox(data, row_start, row_end)
% coordinates and body part indices for frame
coords = data.res(row_start:row_end,2:4) ;
IDX = data.RESIDX(row_start:row_end,:) ;
% get coord rows corresponding to left and right wings
wingRows_R = (IDX(:,data.rightWingInd)==1) ;
wingRows_L = (IDX(:,data.leftWingInd)==1) ;
% get wing coords (in voxels, so int16)
wingVoxR = coords(wingRows_R, :) ;
wingVoxL = coords(wingRows_L, :) ;

end
% -------------------------------------------------------------------------
%OLD CODE TO CHECK MERGED WING CLUSTERING ASSIGNMENT
%{
        dist_mat = pdist2(centroids, [rightWingCM_interp(ind,:) ; ...
            leftWingCM_interp(ind,:)]) ; 
        diagSum = sum(diag(dist_mat)) ; 
        offDiagSum = sum(diag(fliplr(dist_mat))) ; 
        if (offDiagSum - diagSum) > swapThresh
            right_idx = 1 ; 
            left_idx = 2 ; 
            unMergeFailFlag = false ; 
        elseif (offDiagSum - diagSum) < -1*swapThresh
           right_idx = 2 ; 
           left_idx = 1 ; 
           unMergeFailFlag = false ; 
        else
            right_idx = 1 ; 
            left_idx = 1 ;
            unMergeFailFlag = true ; 
        end
        
        min_clust_R_dist = dist_mat(right_idx, 1) ;
        min_clust_L_dist = dist_mat(left_idx, 2) ; 
        % check this against the original error. i.e., if clustering didn't
        % actually change the CM distances much from the original
        % calculation, it might mean clustering was unnecessary in the
        % first place
        right_cm_dist = myNorm(data_in.rightWingCM(ind,:) - ...
            rightWingCM_interp(ind,:)) ;
        left_cm_dist = myNorm(data_in.leftWingCM(ind,:) - ...
            leftWingCM_interp(ind,:)) ;
        new_dist_diff = (right_cm_dist + left_cm_dist) - ...
            (min_clust_R_dist + min_clust_L_dist) ;
        new_dist_check = (new_dist_diff >= distThreshCM/2.0) | ...
            isnan(new_dist_diff) ;
        %}