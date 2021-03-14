%--------------------------------------------------------------------------
% wrapper script for cleanUpWingVoxels
%--------------------------------------------------------------------------
%% path info
% dataPath = 'D:\Fly Data\VNC MN Chrimson\28_05042019\Analysis\Unsorted\Expr_28_mov_003\' ; 
% dataPrefix = 'Expr_28_mov_003' ; 
%dataPath = 'D:\Fly Data\VNC MN Chrimson\03_06112017\Analysis\Unsorted\Expr_3_mov_002\' ; 
%dataPrefix = 'Expr_3_mov_002' ; 
dataPath = 'D:\Fly Data\VNC MN Chrimson\13_25072018\Analysis\Unsorted\Expr_13_mov_014\' ;
dataPrefix = 'Expr_13_mov_014' ; 
suffixStr = '_cleaned' ; 

debugFlag = false ; 
saveFlag = true ; % * note that this only affects saving kinematic plots--data is automatically saved
largePertFlag = false ; 
% ----------------------
%% load important bits
data_in = importdata(fullfile(dataPath, [dataPrefix '_test.mat'])) ; 
% analysisOutput = importdata(fullfile(dataPath, [dataPrefix '_results.mat'])) ; 
data_out_filename = fullfile(dataPath, [dataPrefix suffixStr '.mat']) ; 

% all_fly_bw_xy = analysisOutput.all_fly_bw_xy ; 
% dlt = analysisOutput.dlt_matrix ; 
% dlt_xy = dlt(:, 3) ;                                       
% --------------------------------
%% run function to clean voxels
data = cleanUpWingVoxels(data_in, debugFlag) ; 

% -------------------------------------------
%% check to see if wings need to be swapped
[data, swapFlag] = checkWingSwap(data, largePertFlag) ; 

% ---------------------------
%% re-calculate angles
[rhoTimes, rollVectors] = estimateRollVector(data) ;
data.rhoTimes = rhoTimes ;
data.rollVectors = rollVectors ;

save(data_out_filename, 'data')
data = calcAnglesMain(data_out_filename, largePertFlag, saveFlag, true,...
    suffixStr) ;
