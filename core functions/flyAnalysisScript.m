%==========================================================================
% Script to run analysis taking raw cine files to hull
% reconstruction + position/orientation estimates. Called by
% "flyAnalysisMain.m" and based on "beginning*" scripts
%
% NB: this version of the analysis uses the "batch" tool to run scripts in
% parallel. For sending entire analysis jobs to a cluster, use
% "flyAnalysisScriptClust.m"
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
% allBG   = zeros(3, 512, 512,'uint8') ; % cell(1,3) ;
allBGcell = cell(1,3) ;
allXcm = cell(1,3) ;
allYcm = cell(1,3) ;

% estimations for the time the fly comes in and out of the FOV of each
% camera. this part of the automation can be improved. check it or just set
% the "tin" and "tout" manually later.
allTin  = zeros(3,1) ;
allTout = zeros(3,1) ;

%movieNum = cinFilenames{1}(length(cinFilenames{1})-6:length(cinFilenames{1})-4) ;
movieNum = movNumStr ; 

tic

disp('Finding background images. Batch mode...') ;

job1 = batch('findBG_MOG_mk3',5,{cinFilenames{1}}) ;
job2 = batch('findBG_MOG_mk3',5,{cinFilenames{2}}) ;
job3 = batch('findBG_MOG_mk3',5,{cinFilenames{3}}) ;

wait(job1) ;
wait(job2) ;
wait(job3) ;

%-----------------------------
% error handling
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

out1=fetchOutputs(job1) ; out2=fetchOutputs(job2) ; out3=fetchOutputs(job3) ;
delete(job1) ; delete(job2) ; delete(job3) ;
% allBG(1,:,:) = out1{1} ; allBG(2,:,:) = out2{1} ; allBG(3,:,:) = out3{1} ;
allBGcell{1} = out1{1} ; allBGcell{2} = out2{1} ; allBGcell{3} = out3{1} ;

allTin(1) = out1{2} ; allTin(2) = out2{2} ; allTin(3) = out3{2} ;
allTout(1) = out1{3} ; allTout(2) = out2{3} ; allTout(3) = out3{3} ;

allXcm{1} = out1{4} ; allXcm{2} = out2{4} ; allXcm{3} = out3{4} ;
allYcm{1} = out1{5} ; allYcm{2} = out2{5} ; allYcm{3} = out3{5} ;

clear out1 out2 out3
batchtime = toc ;

disp(['Done finding background images for movie ' movieNum]) ;
disp(batchtime) ;

if isnan(tin) 
    tin = min([max(allTin), 1]) ; % things get funky when the start is after 0
    tout = min(allTout) ;
end
tin = -100;
%-----------------------------
% test background?
%{ 
figure ; 
for file_num = 1:3
    subplot(1,3,file_num)
    imshow(allBGcell{file_num})
    [~, fn_temp, ~] = fileparts(cinFilenames{file_num}) ;
    title(fn_temp,'Interpreter','none')
end
%}
%---------------------------------------------
% save these results in case of error later?
if savePointFlag
   bg_savename = fullfile(savePath, prefixStr, 'BG.mat') ;
   save(bg_savename, 'allBGcell','allXcm','allYcm','allTin',...
       'allTout','tin','tout')
end

%  -----------------------------------------------------------------------
%% PERFORM BINARY THRESHOLDING ON IMAGES
%  -----------------------------------------------------------------------
tic

disp('Doing binary thresholds...')
cam = XY ;
jobXY = batch('binaryThreshold',8,{allBGcell{cam}, ...
    cinFilenames{cam}, tin, tout, twoFlies_xy, allXcm{cam}, allYcm{cam}, ...
    removeLegsFlag, stopWingsFlag}) ;

cam = XZ ;
jobXZ = batch('binaryThreshold',8,{allBGcell{cam}, ...
    cinFilenames{cam}, tin, tout, twoFlies_xz, allXcm{cam}, allYcm{cam}, ...
    removeLegsFlag, stopWingsFlag}) ;

cam = YZ ;
jobYZ = batch('binaryThreshold',8,{allBGcell{cam}, ...
    cinFilenames{cam}, tin, tout, twoFlies_yz, allXcm{cam}, allYcm{cam}, ...
    removeLegsFlag, stopWingsFlag}) ;

% wait on jobs
wait(jobXY) ; wait(jobXZ) ;  wait(jobYZ) ;

%------------------------------
% error handling
if ~isempty(jobXY.Tasks.Error) || ~isempty(jobXZ.Tasks.Error) || ...
        ~isempty(jobYZ.Tasks.Error)
    cam_names_temp = {'xy', 'xz', 'yz'} ;
    jobTasks = [jobXY.Tasks jobXZ.Tasks jobYZ.Tasks] ;
    errors_idx = [~isempty(jobTasks(1).Error) ~isempty(jobTasks(2).Error)...
        ~isempty(jobTasks(3).Error)] ;
    errors_ind = find(errors_idx) ; 
    errors = jobTasks(errors_ind(1)).Error ;
    %disp(getReport(errors, 'basic'))
    msg = strcat('Error doing binary threshold for: ', ...
        [ ' ' cam_names_temp{errors_ind(1)} '_'  movieNum]) ;
    disp(msg)
    msg = strcat(msg, ': ', getReport(errors, 'extended')) ;
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
with_legs_bw_xy       = outXY{8} ; 

all_fly_bw_xz         = outXZ{1} ;
body_only_bw_xz       = outXZ{2} ;
all_fly_thresholds_xz = outXZ{3} ;
xcm_xz                = outXZ{4} ;
ycm_xz                = outXZ{5} ;
allAxlim_xz           = outXZ{6} ;
with_legs_bw_xz       = outXZ{8} ;

all_fly_bw_yz         = outYZ{1} ;
body_only_bw_yz       = outYZ{2} ;
all_fly_thresholds_yz = outYZ{3} ;
xcm_yz                = outYZ{4} ;
ycm_yz                = outYZ{5} ;
allAxlim_yz           = outYZ{6} ;
with_legs_bw_yz       = outYZ{8} ;

clear outXY outXZ outYZ 

disp(['Done with binaryThreshold for movie ' movieNum]) ;
delete(jobXY) ;
delete(jobXZ) ;
delete(jobYZ) ;

binThreshTime = toc ; 
disp(binThreshTime)
%  -----------------------------------------------------------------------
%% COMBINE all_fly_bw_** INTO ONE STRUCTURE
%   (use frames DELTA+1 until Nimages-DELTA)
%  -----------------------------------------------------------------------

dim = max([all_fly_bw_xy.dim ; all_fly_bw_xz.dim ; all_fly_bw_yz.dim]) ; 
%dim = all_fly_bw_xy.dim ;
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
    % ------------------
    % combine XY
    i1=XY ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    % resize XY if need be
    bw_xy = getImage4D(all_fly_bw_xy, 1, i2+DELTA) ;
    if size(bw_xy,1) ~= dim(3) 
        pad_height = round((dim(3) - size(bw_xy,1))/2) ;
        bw_xy = padarray(bw_xy,[pad_height,0],0,'both') ; 
    end
    if size(bw_xy,2) ~= dim(4) 
        pad_width = round((dim(4) - size(bw_xy,2))/2) ;
        bw_xy = padarray(bw_xy,[0, pad_width],0,'both') ; 
    end
    all_fly_bw.mat(ind1vec, ind2vec) = bw_xy ;   
    
    % ------------------
    % combine XZ
    i1=XZ ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    % resize XZ if need be
    bw_xz =  getImage4D(all_fly_bw_xz, 1, i2+DELTA) ; 
    if size(bw_xz,1) ~= dim(3) 
        pad_height = round((dim(3) - size(bw_xz,1))/2) ;
        bw_xz = padarray(bw_xz,[pad_height,0],0,'both') ; 
    end
    if size(bw_xz,2) ~= dim(4) 
        pad_width = round((dim(4) - size(bw_xz,2))/2) ;
        bw_xz = padarray(bw_xz,[0, pad_width],0,'both') ; 
    end
    all_fly_bw.mat(ind1vec, ind2vec) = bw_xz ;
    
    % ------------------
    % combine YZ
    i1=YZ ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    % resize YZ if need be
    bw_yz =   getImage4D(all_fly_bw_yz, 1, i2+DELTA) ;   
    if size(bw_yz,1) ~= dim(3) 
        pad_height = round((dim(3) - size(bw_yz,1))/2) ;
        bw_yz = padarray(bw_yz,[pad_height,0],0,'both') ; 
    end
    if size(bw_yz,2) ~= dim(4) 
        pad_width = round((dim(4) - size(bw_yz,2))/2) ;
        bw_yz = padarray(bw_yz,[0, pad_width],0,'both') ; 
    end
    all_fly_bw.mat(ind1vec, ind2vec) = bw_yz ;
end

% ------------------------------------------------
% combine body_only_bw_** into one structure

dim = max([body_only_bw_xy.dim ; body_only_bw_xz.dim ; ...
    body_only_bw_yz.dim]) ; 
% dim = body_only_bw_xy.dim ;
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
    % --------------------
    % combine XY
    i1=XY ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    % resize XY if need be
    bw_xy = getImage4D(body_only_bw_xy, 1, i2+DELTA) ;
    if size(bw_xy,1) ~= dim(3) 
        pad_height = round((dim(3) - size(bw_xy,1))/2) ;
        bw_xy = padarray(bw_xy,[pad_height,0],0,'both') ; 
    end
    if size(bw_xy,2) ~= dim(4) 
        pad_width = round((dim(4) - size(bw_xy,2))/2) ;
        bw_xy = padarray(bw_xy,[0, pad_width],0,'both') ; 
    end   
    body_only_bw.mat(ind1vec, ind2vec) = bw_xy ;   
    
    % --------------------
    % combine XZ
    i1=XZ ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    % resize XZ if need be
    bw_xz = getImage4D(body_only_bw_xz, 1, i2+DELTA) ; 
    if size(bw_xz,1) ~= dim(3) 
        pad_height = round((dim(3) - size(bw_xz,1))/2) ;
        bw_xz = padarray(bw_xz,[pad_height,0],0,'both') ; 
    end
    if size(bw_xz,2) ~= dim(4) 
        pad_width = round((dim(4) - size(bw_xz,2))/2) ;
        bw_xz = padarray(bw_xz,[0, pad_width],0,'both') ; 
    end
    body_only_bw.mat(ind1vec, ind2vec) = bw_xz;   
    
    % --------------------
    % combine YZ
    i1=YZ ; i2=k ;
    ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
    ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
    % resize YZ if need be
    bw_yz = getImage4D(body_only_bw_yz, 1, i2+DELTA)  ;   
    if size(bw_yz,1) ~= dim(3) 
        pad_height = round((dim(3) - size(bw_yz,1))/2) ;
        bw_yz = padarray(bw_yz,[pad_height,0],0,'both') ; 
    end
    if size(bw_yz,2) ~= dim(4) 
        pad_width = round((dim(4) - size(bw_yz,2))/2) ;
        bw_yz = padarray(bw_yz,[0, pad_width],0,'both') ; 
    end
    body_only_bw.mat(ind1vec, ind2vec) = bw_yz ;   
end

% ------------------------------------------------------
% combine binarized images that include legs, if using

if isfield(with_legs_bw_xy, 'dim')
    dim = max([with_legs_bw_xy.dim ; with_legs_bw_xz.dim ; with_legs_bw_yz.dim]) ; 
    % dim = with_legs_bw_xy.dim ;
    newdim = dim ;
    % newdim(1) = dim(1) - 2*DELTA ;
    % newdim(2) = 3 ; % three cams
    newdim(2) = dim(2) - 2*DELTA ;
    newdim(1) = 3 ; % three cams
    
    with_legs_bw = init4D(newdim) ;
    % Nimages = newdim(1) ;
    Nimages = newdim(2) ;
    for k=1:Nimages
        % -------------------------
        % combine XY
        i1=XY ; i2=k ;
        ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
        ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
        % resize XY if need be
        bw_xy = getImage4D(with_legs_bw_xy, 1, i2+DELTA) ;
        if size(bw_xy,1) ~= dim(3)
            pad_height = round((dim(3) - size(bw_xy,1))/2) ;
            bw_xy = padarray(bw_xy,[pad_height,0],0,'both') ;
        end
        if size(bw_xy,2) ~= dim(4)
            pad_width = round((dim(4) - size(bw_xy,2))/2) ;
            bw_xy = padarray(bw_xy,[0, pad_width],0,'both') ;
        end
        with_legs_bw.mat(ind1vec, ind2vec) = bw_xy ;
        
        % -------------------------
        % combine XZ
        i1=XZ ; i2=k ;
        ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
        ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
        % resize XZ if need be
        bw_xz = getImage4D(with_legs_bw_xz, 1, i2+DELTA) ;
        if size(bw_xz,1) ~= dim(3)
            pad_height = round((dim(3) - size(bw_xz,1))/2) ;
            bw_xz = padarray(bw_xz,[pad_height,0],0,'both') ;
        end
        if size(bw_xz,2) ~= dim(4)
            pad_width = round((dim(4) - size(bw_xz,2))/2) ;
            bw_xz = padarray(bw_xz,[0, pad_width],0,'both') ;
        end
        with_legs_bw.mat(ind1vec, ind2vec) = bw_xz ;
        
        % -------------------------
        % combine YZ
        i1=YZ ; i2=k ;
        ind1vec = dim(3)*(i1-1) +  (1:dim(3)) ;
        ind2vec = dim(4)*(i2-1) +  (1:dim(4)) ;
        % resize YZ if need be
        bw_yz = getImage4D(with_legs_bw_yz, 1, i2+DELTA) ;
        if size(bw_yz,1) ~= dim(3)
            pad_height = round((dim(3) - size(bw_yz,1))/2) ;
            bw_yz = padarray(bw_yz,[pad_height,0],0,'both') ;
        end
        if size(bw_yz,2) ~= dim(4)
            pad_width = round((dim(4) - size(bw_yz,2))/2) ;
            bw_yz = padarray(bw_yz,[0, pad_width],0,'both') ;
        end
        with_legs_bw.mat(ind1vec, ind2vec) = bw_yz ;
    end
else
    with_legs_bw = [] ;
end

disp(['Done combining for movie ' movieNum])
%--------------------------------------------------------------------------
%% COMBINE CENTER-OF-MASS COORDINATES FOR EACH FRAME/CAMERA
%--------------------------------------------------------------------------
% CM_pos is the center-of-mass of the body in each image. 
% dimension of CM_pos is (3cameras, Nimages, 2coordinates)

% first adjust for any size changes we may have made to images
pad_amt_yz = all_fly_bw.dim(3:4) - all_fly_bw_yz.dim(3:4) ; 
pad_amt_xz = all_fly_bw.dim(3:4) - all_fly_bw_xz.dim(3:4) ; 
pad_amt_xy = all_fly_bw.dim(3:4) - all_fly_bw_xy.dim(3:4) ; 

% then combine CM measured from each camera into CM_pos
CM_pos = zeros(3, Nimages, 2) ;
ind = (1:Nimages) + DELTA ;
CM_pos(XY, :,1) = xcm_xy(ind) + pad_amt_xy(2)/2 ;
CM_pos(XY, :,2) = ycm_xy(ind) + pad_amt_xy(1)/2;
CM_pos(XZ, :,1) = xcm_xz(ind) + pad_amt_xz(2)/2;
CM_pos(XZ, :,2) = ycm_xz(ind) + pad_amt_xz(1)/2;
CM_pos(YZ, :,1) = xcm_yz(ind) + pad_amt_yz(2)/2;
CM_pos(YZ, :,2) = ycm_yz(ind) + pad_amt_yz(1)/2;


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
%params.detectorLengthPix=512;
params.detectorLengthPix = all_fly_bw.dim(3:4) ;  % [imageHeight, imageWidth]
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
    'all_fly_thresholds_yz', 'xcm_yz', 'ycm_yz', 'allAxlim_yz', ...
    'with_legs_bw_xy', 'with_legs_bw_xz', 'with_legs_bw_yz', 'all_fly_bw',...
    'body_only_bw', 'with_legs_bw', 'params', 'easyWandData', 'CM_pos')
end

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
%--------------------------------------------------------------------------
%% NEW 3D RECONSTRUCTION
%--------------------------------------------------------------------------
tic ;
disp('Doing hull reconstruction...')
try
    [ bodyRes, bodyFrameStartInd, bodyFrameEndInd, ...
        wing1Res, wing1FrameStartInd, wing1FrameEndInd, ...
        wing2Res, wing2FrameStartInd, wing2FrameEndInd, mergedWingsFlag ] = ...
        hullReconstruction_mk6(params, CM_pos, all_fly_bw, body_only_bw, ...
         dlt_matrix, easyWandData,[2,1,3]);
catch exception
    msg = strcat('Error reconstructing hulls for movie ', movieNum) ;
    msg = strcat(msg, ': ', getReport(exception, 'basic')) ;
    disp(msg)
    fileID = fopen(errorPath,'a+') ;
    fprintf(fileID, '%s\r\n', msg) ;
    fprintf(fileID, '%s\r\n', ' ') ;
    fclose(fileID) ;
    errorFlag = true ;
    delete(gcp) ;
    return
end
% [ bodyRes, bodyFrameStartInd, bodyFrameEndInd, ...
%     wing1Res, wing1FrameStartInd, wingFrameEndInd, ...
%     wing2Res, wing2FrameStartInd, wing2FrameEndInd, mergedWingsFlag ] = ...
%     hullReconstruction_mk6(params, CM_pos, all_fly_bw, body_only_bw, dlt_matrix, easyWandData,[2,1,3]);
thull = toc ;

delete(gcp) ;

disp(['done calculating Hulls for movie ' movieNum]) ;
disp(toc)
%---------------------------------------------
% save these results in case of error later?
if savePointFlag
   hullRecon_savename = fullfile(savePath, prefixStr, 'hullRecon.mat') ;
   save(hullRecon_savename, 'bodyRes', 'bodyFrameStartInd', ...
       'bodyFrameEndInd', 'wing1Res', 'wing1FrameStartInd',...
       'wing1FrameEndInd', 'wing2Res', 'wing2FrameStartInd',...
       'wing2FrameEndInd','mergedWingsFlag','params')
end
%--------------------------------------------------------------------------
%% ANALYZE VOXEL RECONSTRUCTION
%--------------------------------------------------------------------------
t1 = clock ;

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
    tic
    data = hullAnalysis_mk3 (bodyRes, wing1Res, wing2Res, params, ...
        mergedWingsFlag, [], 'test', plotHullFlag, saveHullFigFlag,...
        hullFigPath);
    toc
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
disp(dt12)
diary off ;

%--------------------------------------------------------------------------
%% SAVE RESULTS
%--------------------------------------------------------------------------
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

savePathFull = fullfile(savePath, prefixStr, [resultsFileName '.mat']) ; 
if exist('allBG','var')
save(savePathFull, 'data', 'bodyRes', 'bodyFrameStartInd', 'bodyFrameEndInd', ...
    'wing1Res', 'wing1FrameStartInd', 'wing1FrameEndInd', ...
    'wing2Res', 'wing2FrameStartInd', 'wing2FrameEndInd', 'mergedWingsFlag', ...
    'params',  'CM_pos', 'all_fly_bw', 'body_only_bw',  'with_legs_bw',...
    'dlt_matrix', 'easyWandData', ...
    'cinFilenames', 'allBG', 'Nimages', 'all_fly_bw_xy',  'body_only_bw_xy',...
    'all_fly_thresholds_xy',  'xcm_xy', 'ycm_xy', 'allAxlim_xy', 'DELTA', ...
    'all_fly_bw_xz', 'body_only_bw_xz', 'all_fly_thresholds_xz', 'xcm_xz',...
    'ycm_xz', 'allAxlim_xz', 'all_fly_bw_yz', 'body_only_bw_yz', ...
    'all_fly_thresholds_yz', 'xcm_yz', 'ycm_yz', 'allAxlim_yz', ...
    'tin', 'tout', 'resultsFileName')   
else
save(savePathFull, 'data', 'bodyRes', 'bodyFrameStartInd', 'bodyFrameEndInd', ...
    'wing1Res', 'wing1FrameStartInd', 'wing1FrameEndInd', ...
    'wing2Res', 'wing2FrameStartInd', 'wing2FrameEndInd', 'mergedWingsFlag', ...
    'params',  'CM_pos', 'all_fly_bw', 'body_only_bw',  'with_legs_bw',...
    'dlt_matrix', 'easyWandData', ...
    'cinFilenames', 'allBGcell', 'Nimages', 'all_fly_bw_xy',  'body_only_bw_xy',...
    'all_fly_thresholds_xy',  'xcm_xy', 'ycm_xy', 'allAxlim_xy', 'DELTA', ...
    'all_fly_bw_xz', 'body_only_bw_xz', 'all_fly_thresholds_xz', 'xcm_xz',...
    'ycm_xz', 'allAxlim_xz', 'all_fly_bw_yz', 'body_only_bw_yz', ...
    'all_fly_thresholds_yz', 'xcm_yz', 'ycm_yz', 'allAxlim_yz', ...
    'tin', 'tout', 'resultsFileName')       
end


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
