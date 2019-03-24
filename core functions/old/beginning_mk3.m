%{
phantomSDK_setPath ;
LoadPhantomLibraries ;
%}

%mk2: Trying to make compatible with automation code
%mk3: Speeding up analysis by doing things in batch mode

%% CAMERA CALIBRATION 
%  -----------------------------------------------------------------------

% find points for camera wand calibration
% use:
% parpool(4) ;
% M = getWandPoints (radius1, radius2, xyFilename, xzFilename, yzFilename, outputCSVfile)
% e.g.
%


%{
xyFilenameCalib = [dataPath 'calibration\xy_wand.cin'] ;
xzFilenameCalib = [dataPath 'calibration\xz_wand.cin'] ;
yzFilenameCalib = [dataPath 'calibration\yz_wand.cin'] ;

cinFilenamesCalib = cell(1,3) ;
cinFilenamesCalib{XY} = xyFilenameCalib ;
cinFilenamesCalib{XZ} = xzFilenameCalib ;
cinFilenamesCalib{YZ} = yzFilenameCalib ;

tic
points = getWandPoints([], [],  xyFilenameCalib, xzFilenameCalib, yzFilenameCalib, [dataPath 'calibration\PitchExprCalib.csv']);
toc

% use easyWand to get calibration
% run: easyWand5
% load points csv file
% optim. mode: FL & principal point
% wand span: 0.0195 meters
% no distortions
% initial guess for cameras: 512x512 frame,  f=6000pix, center=(256,256)
% then save calibration results (both the dlt coeffs (11 per cam) and the
% mat file with the entire calibraion info).
%}


%% LOAD CALIBRATION DATA
% dlt_matrix = load(DLT_matrix_CSV_filename) ; % CSV
% load(easyWandData_filename); % contains easyWandData




%% FIND BACKGROUND ETC.
%  -----------------------------------------------------------------------

%[bg, tin, tout, xcm, ycm] = findBG(cinFilename, thresh)

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

tic

disp('Finding background images. Batch mode...') ;

job1 = batch('findBG_MOG_mk3',5,{cinFilenames{1}}) ;
job2 = batch('findBG_MOG_mk3',5,{cinFilenames{2}}) ;
job3 = batch('findBG_MOG_mk3',5,{cinFilenames{3}}) ;

wait(job1) ;
wait(job2) ;
wait(job3) ;

if ~isempty(job1.Tasks.Error) || ~isempty(job2.Tasks.Error) || ~isempty(job3.Tasks.Error)
    jobTasks = [job1.Tasks job2.Tasks job3.Tasks] ;
    errors = [~isempty(jobTasks(1).Error) ~isempty(jobTasks(2).Error) ~isempty(jobTasks(3).Error)] ;
    errors = find(errors) ; errors = jobTasks(errors(1)).Error ;
    %disp(getReport(errors, 'basic'))
    msg = strcat('Error finding background for movie ', movieNum) ;
    disp(msg)
    msg = strcat(msg, ': ', getReport(errors, 'basic')) ;
    fileID = fopen(errorPath,'a+') ;
    fprintf(fileID, '%s\r\n', msg) ;
    fprintf(fileID, '%s\r\n', ' ') ;
    fclose(fileID) ;
    if strcmp(errors.identifier, 'Component:PossibleFalseTrigger') == 1
        falseTriggerFlag = true ;
    end
    errorFlag = true ;
    return
end

out1 = fetchOutputs(job1) ;
out2 = fetchOutputs(job2) ;
out3 = fetchOutputs(job3) ;

delete(job1) ;
delete(job2) ;
delete(job3) ;

allBG(1,:,:) = out1{1} ;
allBG(2,:,:) = out2{1} ;
allBG(3,:,:) = out3{1} ;

allBGcell{1} = out1{1} ;
allBGcell{2} = out2{1} ;
allBGcell{3} = out3{1} ;

allTin(1) = out1{2} ;
allTin(2) = out2{2} ;
allTin(3) = out3{2} ;

allTout(1) = out1{3} ; 
allTout(2) = out2{3} ;
allTout(3) = out3{3} ;

allXcm{1} = out1{4} ;
allXcm{2} = out2{4} ;
allXcm{3} = out3{4} ;

allYcm{1} = out1{5} ;
allYcm{2} = out2{5} ;
allYcm{3} = out3{5} ;

clear out1 out2 out3
batchtime = toc ;

disp(['Done finding background images for movie ' movieNum]) ;
disp(' ') ;

% test background
%{ 
figure ; 
for file_num = 1:3
    subplot(1,3,file_num)
    imshow(allBGcell{file_num})
    [~, fn_temp, ~] = fileparts(cinFilenames{file_num}) ;
    title(fn_temp,'Interpreter','none')
end
%}


%{
tic;
disp('Finding background images...')

for cam = 1:3 
    try
        [bg, tin_curr, tout_curr, xcm_curr, ycm_curr] = findBG_mk5(cinFilenames{cam}); % finds background here.
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
%}
%% 
%DELTA = 24 ; %pretty arbitrary
if isnan(tin) 
    tin = min([max(allTin), 1]) ; % things get funky when the start is after 0
    tout = min(allTout) ;
end

tic

%disp('Running binaryThreshold. batch mode.') ;

%{
figure('Position', [85 400 1650 578]) ; 
subplot(1,3,1) 
imshow(allBGcell{1})
title('YZ')

subplot(1,3,2) 
imshow(allBGcell{2})
title('XZ')

subplot(1,3,3) 
imshow(allBGcell{3})
title('XY')

%}
disp('Doing binary thresholds...')
cam = XY ;
jobXY = batch('binaryThreshold',7,{squeeze(allBG(cam,:,:)), cinFilenames{cam}, tin, tout, twoFlies_xy, allXcm{cam}, allYcm{cam}}) ;
%{
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
%} 

cam = XZ ;
jobXZ = batch('binaryThreshold',7,{squeeze(allBG(cam,:,:)), cinFilenames{cam}, tin, tout, twoFlies_xz, allXcm{cam}, allYcm{cam}}) ;
%{
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
%}

cam = YZ ;
jobYZ = batch('binaryThreshold',7,{squeeze(allBG(cam,:,:)), cinFilenames{cam}, tin, tout, twoFlies_yz, allXcm{cam}, allYcm{cam}}) ;
%{
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
%}
wait(jobXY) ;
wait(jobXZ) ;
wait(jobYZ) ;

if ~isempty(jobXY.Tasks.Error) || ~isempty(jobXZ.Tasks.Error) || ~isempty(jobYZ.Tasks.Error)
    jobTasks = [jobXY.Tasks jobXZ.Tasks jobYZ.Tasks] ;
    errors = [~isempty(jobTasks(1).Error) ~isempty(jobTasks(2).Error) ~isempty(jobTasks(3).Error)] ;
    errors = find(errors) ; errors = jobTasks(errors(1)).Error ;
    %disp(getReport(errors, 'basic'))
    msg = strcat('Error doing binary threshold for movie ', movieNum) ;
    disp(msg)
    msg = strcat(msg, ': ', getReport(errors, 'basic')) ;
    fileID = fopen(errorPath,'a+') ;
    fprintf(fileID, '%s\r\n', msg) ;
    fprintf(fileID, '%s\r\n', ' ') ;
    fclose(fileID) ;
    errorFlag = true ;
    return
end

% assign output arguments and delete jobs
outXY = fetchOutputs(jobXY) ;
outXZ = fetchOutputs(jobXZ) ;
outYZ = fetchOutputs(jobYZ) ;

all_fly_bw_xy         = outXY{1} ;
body_only_bw_xy       = outXY{2} ;
all_fly_thresholds_xy = outXY{3} ;
xcm_xy                = outXY{4} ;
ycm_xy                = outXY{5} ;
allAxlim_xy           = outXY{6} ;
DELTA                 = outXY{7} ;

all_fly_bw_xz         = outXZ{1} ;
body_only_bw_xz       = outXZ{2} ;
all_fly_thresholds_xz = outXZ{3} ;
xcm_xz                = outXZ{4} ;
ycm_xz                = outXZ{5} ;
allAxlim_xz           = outXZ{6} ;

all_fly_bw_yz         = outYZ{1} ;
body_only_bw_yz       = outYZ{2} ;
all_fly_thresholds_yz = outYZ{3} ;
xcm_yz                = outYZ{4} ;
ycm_yz                = outYZ{5} ;
allAxlim_yz           = outYZ{6} ;

clear outXY outXZ outYZ 

disp(['Done with binaryThreshold for movie ' movieNum]) ;
delete(jobXY) ;
delete(jobXZ) ;
delete(jobYZ) ;


toc




%% combine all_fly_bw_** into one structure - used frames DELTA+1 until Nimages-DELTA

dim = all_fly_bw_xy.dim ;
newdim = dim ;
newdim(1) = dim(1) - 2*DELTA ;
newdim(2) = 3 ; % three cams

all_fly_bw = init4D(newdim) ; 
Nimages = newdim(1) ;
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
newdim(1) = dim(1) - 2*DELTA ;
newdim(2) = 3 ; % three cams
body_only_bw = init4D(newdim) ; 
Nimages = newdim(1) ;

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
%%

% combine center-of-mass coordinates for each image in each camera
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


%% VIZ body vs. whole fly featuring
% ---------------------------------
%{
figure('position',[ 94   584   560   160]);

ww = [-1 1 -1 1] * 48  ;

for k=1:Nimages
    for cam=1:3
        fly = getImage4D(all_fly_bw,cam,k) ;
        bod = getImage4D(body_only_bw, cam,k) ;
        rgb = zeros(512,512,3,'uint8') ;        
        rgb(:,:,1) = uint8(bod)*255 ; 
        rgb(:,:,2) = uint8(fly)*255 ; 
        subplot(1,3,cam) ; 
        imshow(rgb) ;
        xc = squeeze(CM_pos(cam,k,1)) ;
        yc = squeeze(CM_pos(cam,k,2)) ;
        axis([xc xc yc yc]+ww) ;
    end
    title(k) ;
    %saveas(gcf,['.\tmp\featuring_' num2str(k) '.png']) ;
    pause (0.05);
end
%}
%% NEW 3D RECONSTRUCTION




tic ;
disp('Doing hull reconstruction...')
try
    [ bodyRes, bodyFrameStartInd, bodyFrameEndInd, ...
        wing1Res, wing1FrameStartInd, wing1FrameEndInd, ...
        wing2Res, wing2FrameStartInd, wing2FrameEndInd, mergedWingsFlag ] = ...
        hullReconstruction_mk5(params, CM_pos, all_fly_bw, body_only_bw, dlt_matrix, easyWandData,[2,1,3]);
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
%% hull analysis
t1 = clock ;

diaryFile = ['myDiary.txt'] ;
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
        mergedWingsFlag, [], 'test');
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

% removing legs/pin
% to do - measure wing vectors with respect to the vein.

if (exist('resultsFileName','var') ~= 1)
    [filepath,name,ext] = fileparts(cinFilenames{1}) ;
    
    name_split = strsplit(name,'_') ; 
    movNumStr_curr = name_split{2} ; 
    
    filepath_split = split(filepath,'\') ; 
    if isempty(filepath_split{end})
        folder_name = filepath_split{end-1} ; 
    else
        folder_name = filepath_split{end} ;
    end
    exprNum_curr = str2double(folder_name(1:2)) ;
    exprNumStr_curr = num2str(exprNum_curr) ;
    
    prefixStr = ['Expr_' exprNumStr_curr '_mov_' movNumStr_curr ];
    resultsFileName = [prefixStr '_results'] ;
end


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


save results_temp data bodyRes bodyFrameStartInd bodyFrameEndInd ...
    wing1Res wing1FrameStartInd wing1FrameEndInd ...
    wing2Res wing2FrameStartInd wing2FrameEndInd mergedWingsFlag ...
    params  CM_pos  all_fly_bw body_only_bw  dlt_matrix  easyWandData ...
    cinFilenames allBG Nimages ...
    all_fly_bw_xy  body_only_bw_xy  all_fly_thresholds_xy  xcm_xy  ycm_xy  allAxlim_xy ...
    DELTA  ...
    all_fly_bw_xz  body_only_bw_xz all_fly_thresholds_xz xcm_xz ycm_xz allAxlim_xz ...
    all_fly_bw_yz  body_only_bw_yz all_fly_thresholds_yz xcm_yz ycm_yz allAxlim_yz ...
    tin tout resultsFileName   

dos(['move/y results_temp.mat ' resultsFileName '.mat']) ;

return
