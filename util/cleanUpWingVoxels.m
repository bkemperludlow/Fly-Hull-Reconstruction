%--------------------------------------------------------------------------
% function to correct errors in wing voxel reconstruction due to occlusion,
% etc. I'm beginning this to deal with cases in which the wings are large
% blobs because e.g. the side view of one wing is being counted as the side
% view of another.
%
% Note: currently, this will only work if most of the frames are okay and
% it's just e.g. the back of the wing stroke that has problems
%--------------------------------------------------------------------------
function data_out = cleanUpWingVoxels(data_in, debugFlag)
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
solidityThresh = 0.8 ; 

% initialize output as copy of input
data_out = data_in ;
% -------------------------------------------------------------------------
%% load in some relevant data
N_frames = data_in.Nimages ;
rightWingInd = data_in.rightWingInd ;
leftWingInd = data_in.leftWingInd ;
bodyCM = data_in.bodyCM ;
wingLength = data_in.wingLength ;


% get row indices
df = diff(data_in.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data_in.res,1)] ;
clear df ;
% -------------------------------------------------------------------------
%% use wing center of mass jumps to determine putative bad recon. frames
[rightWingCM_interp, good_idx_R, rightWingTips_interp] = ...
    interpolateWingCM(data_in,'right') ;
[leftWingCM_interp, good_idx_L, leftWingTips_interp] = ...
    interpolateWingCM(data_in,'left') ;

rightWingTips_interp = filter_wing_vecs(rightWingTips_interp,3,2,2000) ;
leftWingTips_interp = filter_wing_vecs(leftWingTips_interp,3,2,2000) ;

bad_idx = ~good_idx_R | ~good_idx_L ;
bad_frames = find(bad_idx) ;
N_bad_frames = length(bad_frames) ;
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
% -------------------------------------------------------------------------
%% loop through bad frames and see if we can sort out the issue
if debugFlag
    h_debug_in = figure('PaperPositionMode','auto','Position',figPosBefore) ;
    h_debug_out = figure('PaperPositionMode','auto','Position',figPosAfter) ;
end

clusterRightFlag = nan(N_bad_frames,1) ;
clusterLeftFlag = nan(N_bad_frames,1) ;
clusterErrorFlag = nan(N_bad_frames, 1) ;

for k = 1:N_bad_frames %19
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
    % check if we can use what we know to diagnose the failure mode
    rightBadFlag = ~(good_idx_R(ind)) ;
    leftBadFlag = ~(good_idx_L(ind)) ;
    rightVoxFlag = abs(N_vox_R_zScore(ind)) > voxZscoreThresh ;
    leftVoxFlag = abs(N_vox_L_zScore(ind)) > voxZscoreThresh ;
    
    %----------------------------------------------------------------------
    %% if wings are merged, cluster them
    % If all flags are tripped, we guess that 1) the wings
    % are merged abd 2) there's a crappy little bundle of pixels for one
    % wing off in a corner. in this case, we'll cluster
    if (rightBadFlag || leftBadFlag) && (rightVoxFlag && leftVoxFlag)
        fprintf('Frame %d : wings seemed merged, attempting to cluster... \n',...
            ind )
        if N_vox_R_zScore(ind) > N_vox_L_zScore(ind)
            wing_str = 'right' ;
        elseif N_vox_R_zScore(ind) < N_vox_L_zScore(ind)
            wing_str = 'left' ;
        else
            disp('this should not happen')
            keyboard
        end
        % perform clustering
        [wingVox, wingRows, label_idx, centroids, badClusterFlag] = ...
            clusterWings(data_in, ind, wing_str, clustDebugFlag1) ;
        
        if badClusterFlag
            clusterErrorFlag(k) = true ;
            continue
        end
        
        % give the debug plot a title, if it comes up
        if clustDebugFlag1
            title(['Frame ' num2str(ind) ' merged wing clustering'])
        end
        % if we did successfully cluster, we now need to figure out
        % which is left vs right. we'll do this by checking against the
        % interpolated center of mass
        right_centroid_dist = myNorm(centroids -  ...
            repmat(rightWingCM_interp(ind,:),2,1)) ;
        left_centroid_dist = myNorm(centroids -  ...
            repmat(leftWingCM_interp(ind,:),2,1)) ;
        [min_clust_R_dist, right_idx] = min(right_centroid_dist) ;
        [min_clust_L_dist, left_idx] = min(left_centroid_dist) ;
        
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
        new_dist_check = new_dist_diff < distThreshCM/2.0 ;
        
        % in case both blobs think they should be the same wing
        if (right_idx == left_idx) || new_dist_check
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
    % grab wing voxels again, to include changes
    coords = data_out.res(row_start:row_end,2:4) ; 
    IDX = data_out.RESIDX(row_start:row_end,:) ; 
    wingRows_R = (IDX(:,rightWingInd)==1) ;
    wingRows_L = (IDX(:,leftWingInd)==1) ;
    wingVoxR = coords(wingRows_R, :) ; 
    wingVoxL = coords(wingRows_L, :) ; 
    
    % perform some checks to see if we should cluster
    checkRightCM = (myNorm(data_out.rightWingCM(ind,:) - ...
        rightWingCM_interp(ind,:)) > distThreshCM) ;
    checkLeftCM = (myNorm(data_out.leftWingCM(ind,:) - ...
        leftWingCM_interp(ind,:)) > distThreshCM) ;
    
    justRightCheck = checkRightCM | (N_vox_R_zScore(ind) > voxZscoreThresh) ;
    justLeftCheck = checkLeftCM | (N_vox_L_zScore(ind) > voxZscoreThresh) ;
    
    [rightWingSolidity, ~, ~] = wingVoxSolidity(wingVoxR) ;
    [leftWingSolidity, ~, ~] = wingVoxSolidity(wingVoxL) ;
    
    solidityCheckR = (rightWingSolidity > solidityThresh) & ...
        (abs(N_vox_R_zScore(ind)) < voxZscoreThresh) ;
    solidityCheckL = (leftWingSolidity > solidityThresh) & ...
        (abs(N_vox_L_zScore(ind)) < voxZscoreThresh) ;
    % if both of these are true, either our interpolation is bad, or the
    % clustering wasn't good in the first place--need to add this back in?
    try
        oneWingClustFlag = false ; 
        if justRightCheck && ~solidityCheckR
            fprintf('Frame %d : clustering just right wing voxels \n', ind )
            wing_str = 'right' ;
            [data_out, ~] = clusterSingleWing(data_out, wing_str, ind, ...
                rightWingCM_interp, rightWingTips_interp, row_start, row_end,...
                clustDebugFlag2) ;
            oneWingClustFlag = true ; 
        end
        if justLeftCheck && ~solidityCheckL
            fprintf('Frame %d : clustering just left wing voxels \n', ind )
            wing_str = 'left' ;
            [data_out, ~] = clusterSingleWing(data_out, wing_str, ind, ...
                leftWingCM_interp, leftWingTips_interp, row_start, row_end,...
                clustDebugFlag2) ;
            oneWingClustFlag = true ; 
        end
        if ~oneWingClustFlag
            fprintf('Frame %d : no blobby wings to cluster \n', ind)
        end
        
    catch
        keyboard
    end
    %----------------------------------------------------------------------
    % generate plots to check our work?
    if debugFlag
        h_debug_out = plotFlyVoxels(data_out, ind, h_debug_out) ;
        title(['Frame ' num2str(ind) ' AFTER'])
    end
    disp(k)
end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper functions
function [data_out, badClusterFlag] = clusterSingleWing(data_in, wing_side,...
    ind, wingCM_interp, wingTip_interp, row_start, row_end, clusterDebugFlag)
%--------------------------
% inputs and params
% data_out starts as a copy of data_in
data_out = data_in ;

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
end

%------------------------------------------
% determine which cluster gives best wing
% NB: just taking the longest for now, should fix soon. maybe compare
% against XY view BW image?
tip2span_dists = myNorm(wingTipMat - centroids) ;
cm_dists = myNorm(repmat(wingCM_interp(ind,:),N_clusts,1) - centroids) ;
pca_ratio = cellfun(@(y) y(1)/y(3), pcaLatentCell) ; 
pca_span_dot = cellfun(@(y) abs(dot(spanPrev,y(:,1))), pcaCoeffCell) ;
% pca_axis_cross_z = cellfun(@(y) y(3), pca_axis_cross) ;

[~, max_dist_idx] = max(tip2span_dists) ;
[~, closest_cm_idx] = min(cm_dists) ; 
[~, max_pca_idx] = max(pca_ratio) ;
% if abs(tip2span_dists(1) - tip2span_dists(2)) > tipDistThresh
%     good_idx = max_dist_idx ; 
if (abs(pca_ratio(1) - pca_ratio(2)) > pcaRatioThresh) && ...
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
if (sum(label_idx == good_idx) < minNumVox) || (abs(dotCheck) < dotThresh)
    disp('bad clustering--just take original')
else
    %----------------------------------
    % else enter data into structure
    data_out.([wing_side 'WingCM'])(ind,:) = centroids(good_idx,:) ;
    data_out.([wing_side 'SpanHats'])(ind,:) = spanHatMat(good_idx, :) ;
    data_out.([wing_side 'ChordHats'])(ind,:) = chordHatMat(good_idx,:) ;
    data_out.([wing_side 'ChordAltHats'])(ind,:) = chordAltHatMat(good_idx,:) ;
    data_out.([wing_side 'WingTips'])(ind,:) = wingTipMat(good_idx,:) ;
    data_out.RESIDX(row_start:row_end,wingInd) = wingRows(:, good_idx) ;
end

end