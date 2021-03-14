%--------------------------------------------------------------------------
% function to run 'removeHullGlobs' on the wings with data structure as
% input so we can save the cleaned wing voxels for later use
%--------------------------------------------------------------------------
function data_out = removeWingVoxNoise(data_in, ind_list, verboseFlag, debugFlag)
%% inputs/params
if ~exist('verboseFlag','var')
    verboseFlag = false ;
end
if ~exist('debugFlag','var')
    debugFlag = false ;
end
%% initialize data_out
data_out = data_in ;
rightWingInd = data_in.rightWingInd ;
leftWingInd = data_in.leftWingInd ;

% -----------------------------------------------------------
%% get rows in data structure corresponding to frame "ind"
df = diff(data_in.res(:,1)) ;
frameStartInd = [1 ; find(df==1)+1] ;
frameEndInd   = [frameStartInd(2:end)-1 ; size(data_in.res,1)] ;
clear df ;

for ind = ind_list
    row_start = frameStartInd(ind) ; % rows for voxel idx array
    row_end = frameEndInd(ind) ;
    
    % -------------------------------------
    %% get voxel coordinates for frame
    coords = data_in.res(row_start:row_end,2:4) ;
    IDX = data_in.RESIDX(row_start:row_end,:) ;
    wingRows_R = (IDX(:,rightWingInd)==1) ;
    wingRows_L = (IDX(:,leftWingInd)==1) ;
    wingVoxR = coords(wingRows_R, :) ;
    wingVoxL = coords(wingRows_L, :) ;
    
    % ------------------------
    %% clean wing voxels
    [wingVoxR_cleaned, ~] = removeHullGlobs(wingVoxR) ;
    [wingVoxL_cleaned, ~] = removeHullGlobs(wingVoxL) ;
    
    % -------------------------------------
    %% add cleaned hulls to data structure
    % ------------------------------
    data_out = assignCleanVoxelsToStruct(data_out, row_start, ...
        row_end, wingVoxR, wingVoxR_cleaned, wingRows_R, rightWingInd)  ;
    
    % ------------
    % left wing
    data_out = assignCleanVoxelsToStruct(data_out, row_start, ...
        row_end, wingVoxL, wingVoxL_cleaned, wingRows_L, leftWingInd)  ;
    
    % -------------------------------------
    %% plot results?
    if debugFlag
        h_before = plotFlyVoxels(data_in, ind) ;
        title('before')
        h_after = plotFlyVoxels(data_out, ind) ;
        title('after')
        pause(0.1)
    end
    
    % -------------------------------------
    %% print progress?
    if verboseFlag
        fprintf('Completed %d / %d frames \n',ind, length(ind_list))
    end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% helper function because i'll be damned if i'm writing this shit out twice
function data_out = assignCleanVoxelsToStruct(data_out, row_start, ...
    row_end, wingVox, wingVox_cleaned, wingRows, wingInd)   
% -------------------------------------------------------------------------
% *method taken from:
% https://www.mathworks.com/matlabcentral/answers/196846-finding-matching-rows-in-matrix-for-multiple-points-without-a-loop
    
% ---------
% step 1: combine original and clean voxel arrays and then sort them. since
% we're sorting by rows, identical ones will end up next to each other
% (i.e. the row j and row j+1 will be the same, with j coming from the
% original and j+1 coming from the cleaned voxels)
[vox_sort, vox_sort_idx] = sortrows([wingVox;  wingVox_cleaned]) ;

% ---------
% step 2: find the indices for adjacent matching rows
vox_comb_match_idx = all(vox_sort(1:end-1,:)==vox_sort(2:end,:),2);

% ---------
% step 3: convert the indices from the combined voxel array into indices just
% for the original voxel array
vox_match_sub = vox_sort_idx(vox_comb_match_idx) ;

% ---------
% step 4: because of nested fucking indexing, we're going to find the
% entries in residx that are true for the given wing (since residx has a 
% row for all coordinates). Then we're going to eliminate the ones that
% don't match the cleaned version, and create a new wingRows array
wingRows_sub = find(wingRows) ; % get linear position
wingRows_sub_cleaned = wingRows_sub(vox_match_sub) ; % take only matching voxels (i.e. ones in cleaned and original)
wingRows_cleaned = false(size(wingRows)) ;  % initialize new, cleaned wingRows
wingRows_cleaned(wingRows_sub_cleaned) = true ; % mark cleaned voxels true

% ------------
% step 5: assign to data struct
data_out.RESIDX(row_start:row_end,wingInd) = wingRows_cleaned ; 
end