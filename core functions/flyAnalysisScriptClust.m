%==========================================================================
% Script to run analysis taking raw cine files to hull
% reconstruction + position/orientation estimates. Called by
% "flyAnalysisMain.m" and based on "beginning*" scripts
%
% NB: this version of the analysis runs through things linearly (i.e. not 
% parallel). For running analysis of a movie using the batch function, use
% "flyAnalysisScript.m"
%
% The basic chain of events is:
%   1) background estimation for each camera view (findBG_MOG_mk3.m)
%   2) binarization of images from each camera that tries to separate body
%   vs non-body (wing) pixels (binaryThreshold.m)
%   3) generate voxel reconstruction of fly at each frame based on
%   binarized images and camera calibration (hullReconstruction_mk5.m)
%   4) estimate vectors that describe the fly's DOF based on voxel
%   representation, e.g. body axis vectors, wing span vectors, etc
%   (hullAnalysis_mk3,m)
%==========================================================================

%% FIND BACKGROUND ETC.
%  -----------------------------------------------------------------------
% background images in two formats (matrix and cell). I probably only use
% the matrix form.
allBG   = zeros(3, 512, 512,'uint8') ; % cell(1,3) ;
allBGcell = cell(1,3) ;
allXcm = cell(1,3) ;
allYcm = cell(1,3) ;

% estimations for the time the fly comes in and out of the FOV of each
% camera. this part of the automation can be improved. check it or just set
% the "tin" and "tout" manually later.
allTin  = zeros(3,1) ;
allTout = zeros(3,1) ;

movieNum = cinFilenames{1}(length(cinFilenames{1})-6:length(cinFilenames{1})-4) ;

tic;
disp('Finding background images...')

for cam = 1:3 
    try
        [bg, tin_curr, tout_curr, xcm_curr, ycm_curr] = findBG_MOG_mk3(cinFilenames{cam}); % finds background here.
    catch exception
        msg = strcat('Error finding background for movie ', movieNum) ;%cinFilenames{cam}(length(cinFilenames{cam})-6:length(cinFilenames{cam})-4)) ;
        msg = strcat(msg, ': ', getReport(exception, 'basic')) ;
        disp(msg)
        fileID = fopen(errorPath,'a+') ;
        fprintf(fileID, '%s\r\n', msg) ;
        fclose(fileID) ;
        if strcmp(exception.identifier, 'Component:PossibleFalseTrigger') == 1
            falseTriggerFlag = true ;
        end
        %disp(exception.identifier)
        %disp(cam)
%         while cam < 3
%             temp = myOpenCinFile(cinFilenames{cam+1}) ;
%             myCloseCinFile(temp) ;
%             cam = cam + 1 ;
%         end
        return
    end

    allBG(cam,:,:) = bg ;
    allBGcell{cam} = bg ;
    allXcm{cam} = xcm_curr;
    allYcm{cam} = ycm_curr ;
    allTin(cam) = tin_curr ;
    allTout(cam) = tout_curr ;    
end

regtime = toc ;
disp('Done finding background images...') ;

if isnan(tin) 
    tin = max(allTin) ;
    tout = min(allTout) ;
end
%---------------------------------------------
% save these results in case of error later?
if savePointFlag
   bg_savename = fullfile(savePath, prefixStr, 'BG.mat') ;
   save(bg_savename, 'allBG','allBGcell','allXcm','allYcm','allTin',...
       'allTout','tin','tout')
end
%  -----------------------------------------------------------------------
%% PERFORM BINARY THRESHOLDING ON IMAGES
%  -----------------------------------------------------------------------
tic

disp('Doing binary thresholds...')
cam = XY ;
%jobXY = batch('binaryThreshold',7,{squeeze(allBG(cam,:,:)), cinFilenames{cam}, tin, tout, twoFlies_xy, allXcm{cam}, allYcm{cam}}) ;
try
    [all_fly_bw_xy, body_only_bw_xy, all_fly_thresholds_xy, xcm_xy, ycm_xy, allAxlim_xy, DELTA] = ...
        binaryThreshold( squeeze(allBG(cam,:,:)) , cinFilenames{cam}, tin, tout,twoFlies_xy, allXcm{cam}, allYcm{cam}) ;
catch exception
    msg = strcat('Error doing xy binary threshold for movie ', movieNum) ;%cinFilenames{cam}(length(cinFilenames{cam})-6:length(cinFilenames{cam})-4)) ;
    msg = strcat(msg, ': ', getReport(exception, 'basic')) ;
    disp(msg)
    fileID = fopen(errorPath,'a+') ;
    fprintf(fileID, '%s\r\n', msg) ;
    fclose(fileID) ;
    errorFlag = true ;
    return
end

cam = XZ ;
%jobXZ = batch('binaryThreshold',7,{squeeze(allBG(cam,:,:)), cinFilenames{cam}, tin, tout, twoFlies_xz, allXcm{cam}, allYcm{cam}}) ;
try
    [all_fly_bw_xz, body_only_bw_xz, all_fly_thresholds_xz, xcm_xz, ycm_xz, allAxlim_xz, DELTA] = ...
        binaryThreshold( squeeze(allBG(cam,:,:)) , cinFilenames{cam}, tin, tout,twoFlies_xz, allXcm{cam}, allYcm{cam}) ;
catch exception
    msg = strcat('Error doing xz binary threshold for movie ', movieNum) ;%cinFilenames{cam}(length(cinFilenames{cam})-6:length(cinFilenames{cam})-4)) ;
    msg = strcat(msg, ': ', getReport(exception, 'basic')) ;
    fileID = fopen(errorPath,'a+') ;
    fprintf(fileID, '%s\r\n', msg) ;
    fclose(fileID) ;
    errorflag = true ;
    return
end

cam = YZ ;
%jobYZ = batch('binaryThreshold',7,{squeeze(allBG(cam,:,:)), cinFilenames{cam}, tin, tout, twoFlies_yz, allXcm{cam}, allYcm{cam}}) ;
try
    [all_fly_bw_yz, body_only_bw_yz, all_fly_thresholds_yz, xcm_yz, ycm_yz, allAxlim_yz, DELTA] = ...
        binaryThreshold( squeeze(allBG(cam,:,:)) , cinFilenames{cam}, tin, tout,twoFlies_yz, allXcm{cam}, allYcm{cam}) ;
catch exception
    msg = strcat('Error doing yz binary threshold for movie ', movieNum) ;%cinFilenames{cam}(length(cinFilenames{cam})-6:length(cinFilenames{cam})-4)) ;
    msg = strcat(msg, ': ', getReport(exception, 'basic')) ;
    disp(msg)
    fileID = fopen(errorPath,'a+') ;
    fprintf(fileID, '%s\r\n', msg) ;
    fclose(fileID) ;
    errorFlag = true ;
    return
end

toc

%  -----------------------------------------------------------------------
%% COMBINE all_fly_bw_** INTO ONE STRUCTURE
%   (use frames DELTA+1 until Nimages-DELTA)
%  -----------------------------------------------------------------------

dim = all_fly_bw_xy.dim ;
newdim = dim ;
% newdim(1) = dim(1) - 2*DELTA ;
% newdim(2) = 3 ; % three cams
newdim(2) = dim(2) - 2*DELTA ;
newdim(1) = 3 ; % three cams

all_fly_bw = init4D(newdim) ; 
% Nimages = newdim(1) ;
Nimages = newdim(2) ;
disp('Combining...')
for k=1:Nimages
    % combine XY
    i1=XY ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    all_fly_bw.mat(ind1vec, ind2vec) = getImage4D(all_fly_bw_xy, 1, i2+DELTA) ;   
    
    % combine XZ
    i1=XZ ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    all_fly_bw.mat(ind1vec, ind2vec) = getImage4D(all_fly_bw_xz, 1, i2+DELTA) ;   
    
    % combine YZ
    i1=YZ ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    all_fly_bw.mat(ind1vec, ind2vec) = getImage4D(all_fly_bw_yz, 1, i2+DELTA) ;   
end


% combine body_only_bw_** into one structure

dim = body_only_bw_xy.dim ;
newdim = dim ;
% newdim(1) = dim(1) - 2*DELTA ;
% newdim(2) = 3 ; % three cams
% body_only_bw = init4D(newdim) ; 
% Nimages = newdim(1) ;
newdim(2) = dim(2) - 2*DELTA ;
newdim(1) = 3 ; % three cams
body_only_bw = init4D(newdim) ; 
Nimages = newdim(2) ;

% WHEN dealing with body-only need to handle fucking delta.

for k=1:Nimages
    % combine XY
    i1=XY ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    body_only_bw.mat(ind1vec, ind2vec) = getImage4D(body_only_bw_xy, 1, i2+DELTA) ;   
    
    % combine XZ
    i1=XZ ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    body_only_bw.mat(ind1vec, ind2vec) = getImage4D(body_only_bw_xz, 1, i2+DELTA) ;   
    
    % combine YZ
    i1=YZ ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    body_only_bw.mat(ind1vec, ind2vec) = getImage4D(body_only_bw_yz, 1, i2+DELTA) ;   
end
disp(['Done combining for movie ' movieNum])

%--------------------------------------------------------------------------
%% COMBINE CENTER-OF-MASS COORDINATES FOR EACH FRAME/CAMERA
%--------------------------------------------------------------------------
% CM_pos is the center-of-mass of the body in each image. 
% dimension of CM_pos is (3cameras, Nimages, 2coordinates)

CM_pos = zeros(3, Nimages, 2) ;
ind = (1:Nimages) + DELTA ;
CM_pos(XY, :,1) = xcm_xy(ind) ;
CM_pos(XY, :,2) = ycm_xy(ind) ;
CM_pos(XZ, :,1) = xcm_xz(ind) ;
CM_pos(XZ, :,2) = ycm_xz(ind) ;
CM_pos(YZ, :,1) = xcm_yz(ind) ;
CM_pos(YZ, :,2) = ycm_yz(ind) ;


%% DEFINE PARAMS
% --------------

params.CAMERAS=[1 2 3];
params.NCAMS=3;
params.YZ = YZ ; 
params.XY = XY ;
params.XZ = XZ ;
params.fps = 8000 ;
% params.camerasPos=[easyWandData.DLTtranslationVector(:,:,2)';...
% easyWandData.DLTtranslationVector(:,:,1)';...
% easyWandData.DLTtranslationVector(:,:,3)'];  % row1=yz, row2=xz, row3=xy
params.cameraNames = ['yz';'xz';'xy'];
params.detectorLengthPix=512;
params.voxelSize = 50e-6 ; % 50 microns
%params.N=120; % 
params.volLength= 8e-3 ; % size of the square sub-vol cube to reconstruct (meters)
%params.voxelSize=0.4/120;
params.volCenter=[0,0,0];
params.focusPix=easyWandData.focalLengths;
%params.offsetsMatrix=[];
params.startTrackingTime   =  tin+DELTA  ; % plus 1 removed.
params.endTrackingTime     =  tout-DELTA ;
params.firstTrackableFrame =  tin+DELTA  ; % plus 1 removed.

% this is in fact "voxels per cm".
params.pixPerCM            = 350 ; % 232 ; % effective value. need to change that to be consistent with real units
% probably we should start from the wing legnth in cm and calculate how
% many voxels are in 1 wing length by dividing wingLcm/params.voxelSize

%---------------------------------------------
% save these results in case of error later?
if savePointFlag
   binaryThresh_savename = fullfile(savePath, prefixStr, 'binaryThresh.mat') ;
   save(binaryThresh_savename, 'all_fly_bw_xy',  'body_only_bw_xy',...
    'all_fly_thresholds_xy',  'xcm_xy', 'ycm_xy', 'allAxlim_xy', 'DELTA', ...
    'all_fly_bw_xz', 'body_only_bw_xz', 'all_fly_thresholds_xz', 'xcm_xz',...
    'ycm_xz', 'allAxlim_xz', 'all_fly_bw_yz', 'body_only_bw_yz', ...
    'all_fly_thresholds_yz', 'xcm_yz', 'ycm_yz', 'allAxlim_yz')
end
%% VIZ body vs. whole fly featuring
% ---------------------------------
% figure('position',[ 94   584   560   160]);
% 
% ww = [-1 1 -1 1] * 48  ;
% 
% for k=1:Nimages
%     for cam=1:3
%         fly = getImage4D(all_fly_bw,cam,k) ;
%         bod = getImage4D(body_only_bw, cam,k) ;
%         rgb = zeros(512,512,3,'uint8') ;        
%         rgb(:,:,1) = uint8(bod)*255 ; 
%         rgb(:,:,2) = uint8(fly)*255 ; 
%         subplot(1,3,cam) ; 
%         imshow(rgb) ;
%         xc = squeeze(CM_pos(cam,k,1)) ;
%         yc = squeeze(CM_pos(cam,k,2)) ;
%         axis([xc xc yc yc]+ww) ;
%     end
%     title(k) ;
%     %saveas(gcf,['.\tmp\featuring_' num2str(k) '.png']) ;
%     pause (0.05);
% end

%--------------------------------------------------------------------------
%% NEW 3D RECONSTRUCTION
%--------------------------------------------------------------------------
tic ;
disp('Doing hull reconstruction...')
try
    [ bodyRes, bodyFrameStartInd, bodyFrameEndInd, ...
        wing1Res, wing1FrameStartInd, wing1FrameEndInd, ...
        wing2Res, wing2FrameStartInd, wing2FrameEndInd, mergedWingsFlag ] = ...
        hullReconstruction_mk5(params, CM_pos, all_fly_bw, body_only_bw,...
        dlt_matrix, easyWandData,[2,1,3]);
catch exception
    msg = strcat('Error reconstructing hulls for movie ', movieNum) ;%cinFilenames{cam}(length(cinFilenames{cam})-6:length(cinFilenames{cam})-4)) ;
    msg = strcat(msg, ': ', getReport(exception, 'basic')) ;
    disp(msg)
    fileID = fopen(errorPath,'a+') ;
    fprintf(fileID, '%s\r\n', msg) ;
    fprintf(fileID, '%s\r\n', ' ') ;
    fclose(fileID) ;
    errorFlag = true ;
    return
end
% [ bodyRes, bodyFrameStartInd, bodyFrameEndInd, ...
%     wing1Res, wing1FrameStartInd, wingFrameEndInd, ...
%     wing2Res, wing2FrameStartInd, wing2FrameEndInd, mergedWingsFlag ] = ...
%     hullReconstruction_mk6(params, CM_pos, all_fly_bw, body_only_bw, dlt_matrix, easyWandData,[2,1,3]);
thull = toc ;

delete(gcp) ;

disp(['done calculating Hulls for movie ' movieNum]) ;

%---------------------------------------------
% save these results in case of error later?
if savePointFlag
   hullRecon_savename = fullfile(savePath, prefixStr, 'hullRecon.mat') ;
   save(hullRecon_savename, 'bodyRes', 'bodyFrameStartInd', ...
       'bodyFrameEndInd', 'wing1Res', 'wing1FrameStartInd',...
       'wing1FrameEndInd', 'wing2Res', 'wing2FrameStartInd',...
       'wing2FrameEndInd')
end
%--------------------------------------------------------------------------
%% ANALYZE VOXEL RECONSTRUCTION
%--------------------------------------------------------------------------
t1 = clock ;

%diaryFile = ['myDiary.txt'] ;
diaryFile = fullfile(savePath,prefixStr,'myDiary.txt');
try
    dos(['del ' diaryFile ]) ;
catch
    disp('Diary file does not exist.') ;
    pause(2) ;
end
diary(diaryFile) ;
disp('Doing hull analysis...') 
try
    data = hullAnalysis_mk3 (bodyRes, wing1Res, wing2Res, params, ...
        mergedWingsFlag, [], 'test', plotHullFlag, saveHullFigFlag, ...
        hullFigPath);
catch exception
    msg = strcat('Error analyzing hulls for movie ', movieNum) ;%cinFilenames{cam}(length(cinFilenames{cam})-6:length(cinFilenames{cam})-4)) ;
    msg = strcat(msg, ': ', getReport(exception, 'basic')) ;
    disp(msg)
    fileID = fopen(errorPath,'a+') ;
    fprintf(fileID, '%s\r\n', msg) ;
    fprintf(fileID, '%s\r\n', ' ') ;
    fclose(fileID) ;
    errorFlag = true ;
    return
end
disp(['Done with hull analysis for movie ' movieNum])
t2 = clock ;
dt12 = t2 - t1 ;

diary off ;
%--------------------------------------------------------------------------
%% SAVE RESULTS
%--------------------------------------------------------------------------

savePathFull = fullfile(savePath, prefixStr, [resultsFileName '.mat']) ; 
save(savePathFull, 'data', 'bodyRes', 'bodyFrameStartInd', 'bodyFrameEndInd', ...
    'wing1Res', 'wing1FrameStartInd', 'wing1FrameEndInd', ...
    'wing2Res', 'wing2FrameStartInd', 'wing2FrameEndInd', 'mergedWingsFlag', ...
    'params',  'CM_pos', 'all_fly_bw', 'body_only_bw',  'dlt_matrix', 'easyWandData', ...
    'cinFilenames', 'allBG', 'Nimages', 'all_fly_bw_xy',  'body_only_bw_xy',...
    'all_fly_thresholds_xy',  'xcm_xy', 'ycm_xy', 'allAxlim_xy', 'DELTA', ...
    'all_fly_bw_xz', 'body_only_bw_xz', 'all_fly_thresholds_xz', 'xcm_xz',...
    'ycm_xz', 'allAxlim_xz', 'all_fly_bw_yz', 'body_only_bw_yz', ...
    'all_fly_thresholds_yz', 'xcm_yz', 'ycm_yz', 'allAxlim_yz', ...
    'tin', 'tout', 'resultsFileName')   
%dos(['move/y results_temp.mat ' resultsFileName '.mat']) ;

return

%{
% OLD CODE BEFORE hullReconstruction_mk5, maybe you'll need it at some
point (?)

%% PERFORM 3D RECONSRUCTION OF BODY-ONLY
%  -------------------------------------

% there is a problem in finding the body-only image on the xy view. it is
% somehow shorter than it should be.
tic
disp('Find body-only hull') ;
[ bodyHull , bodyHullFrameStartInd, bodyHullFameEndInd] = hullReconstruction_mk4(params, CM_pos, body_only_bw, dlt_matrix, easyWandData);
tbody = toc ;

% PERFORM 3D RECONSRUCTION OF ENTIRE FLY 
%  --------------------------------------
tic
disp('Find entire fly hull') ;
[ flyHull , flyHullFrameStartInd, flyHullFameEndInd] = hullReconstruction_mk4(params, CM_pos, all_fly_bw, dlt_matrix, easyWandData);
tfly = toc ;

%}
