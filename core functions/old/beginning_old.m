
XZ = 2 ;
XY = 3 ;
YZ = 1 ;

%dataPath = 'F:\Pitch\06_220514\' ;
%dataPath = 'F:\Pitch\07_230514\' ;
dataPath = 'C:\Users\Fruit Flies\Desktop\temp\' ;
DLT_matrix_CSV_filename = [dataPath 'calibration_dltCoefs.csv'] ;
easyWandData_filename = [dataPath 'calibration_easyWandData.mat'] ;

format compact


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
dlt_matrix = load(DLT_matrix_CSV_filename) ; % CSV
load(easyWandData_filename); % contains easyWandData


%% DEFINE CIN FILE NAMES

cinFilenames = cell(1,3) ;
cinFilenames{XY} = [dataPath 'xy_009.cin'] ;
cinFilenames{XZ} = [dataPath 'xz_009.cin'] ;
cinFilenames{YZ} = [dataPath 'yz_009.cin'] ;


%% FIND BACKGROUND ETC.
%  -----------------------------------------------------------------------

%[bg, tin, tout, xcm, ycm] = findBG(cinFilename, thresh)

% background images in two formats (matrix and cell). I probably only use
% the matrix form.
allBG   = zeros(3, 512, 512,'uint8') ; % cell(1,3) ;
allBGcell = cell(1,3) ;

% estimations for the time the fly comes in and out of the FOV of each
% camera. this part of the automation can be improved. check it or just set
% the "tin" and "tout" manually later.
allTin  = zeros(3,1) ;
allTout = zeros(3,1) ;


for cam = 1:3 
    [bg, tin, tout] = findBG(cinFilenames{cam}); % finds background here.
    allBG(cam,:,:) = bg ;
    allBGcell{cam} = bg ;
    allTin(cam) = tin ;
    allTout(cam) = tout ;    
end

%% 

% find common tin and tout
tin = -600 ; % max(allTin) ;
tout = 400 ; % min(allTout) ;

cam = XY ;
[all_fly_bw_xy, body_only_bw_xy, all_fly_thresholds_xy, xcm_xy, ycm_xy, allAxlim_xy, DELTA] = ...
    binaryThreshold( squeeze(allBG(cam,:,:)) , cinFilenames{cam}, tin, tout) ;

cam = XZ ;
[all_fly_bw_xz, body_only_bw_xz, all_fly_thresholds_xz, xcm_xz, ycm_xz, allAxlim_xz, DELTA] = ...
    binaryThreshold( squeeze(allBG(cam,:,:)) , cinFilenames{cam}, tin, tout) ;

cam = YZ ;
[all_fly_bw_yz, body_only_bw_yz, all_fly_thresholds_yz, xcm_yz, ycm_yz, allAxlim_yz, DELTA] = ...
    binaryThreshold( squeeze(allBG(cam,:,:)) , cinFilenames{cam}, tin, tout) ;


%% combine all_fly_bw_** into one structure - used frames DELATA+1 until Nimages-DELTA

dim = all_fly_bw_xy.dim ;
newdim = dim ;
newdim(1) = dim(1) - 2*DELTA ;
newdim(2) = 3 ; % three cams

all_fly_bw = init4D(newdim) ; 
Nimages = newdim(1) ;

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

%%

% combine center-of-mass coordinates for each image in each camera
% CM_pos is the center-of-mass of the body in each image. 
% dimension of CM_pos is (3camera, Nimages, 2coordinates)

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
% params.camerasPos=[easyWandData.DLTtranslationVector(:,:,2)';...
% easyWandData.DLTtranslationVector(:,:,1)';...
% easyWandData.DLTtranslationVector(:,:,3)'];  % row1=yz, row2=xz, row3=xy
params.cameraNames = ['yz';'xz';'xy'];
params.detectorLengthPix=512;
params.voxelSize = 50e-6 ; % 50 microns
%params.N=120; % 
params.volLength= 7e-3 ; % size of the square sub-vol cube to reconstruct (meters)
%params.voxelSize=0.4/120;
params.volCenter=[0,0,0];
params.focusPix=easyWandData.focalLengths;
%params.offsetsMatrix=[];
params.startTrackingTime   =  tin+DELTA  ;
params.endTrackingTime     =  tout-DELTA ;
params.firstTrackableFrame =  tin+DELTA  ;

% this is in fact "voxels per cm".
params.pixPerCM            = 232 ; % effective value. need to change that to be consistent with real units
% probably we should start from the wing legnth in cm and calculate how
% many voxels are in 1 wing length by dividing wingLcm/params.voxelSize


%% VIZ body vs. whole fly featuring
% ---------------------------------
figure('position',[ 94   584   560   160]);

ww = [-1 1 -1 1] * 40  ;
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
    %saveas(gcf,['.\tmp\featuring_' num2str(k) '.png']) ;
    pause(0.1) ;
end

%% NEW 3D RECONSTRUCTION

disp('need to add an offset in hullReconstruction_mk5 to translate the reconstructed sub-vol to real-real space') ;

tic ;
[ bodyRes, bodyFrameStartInd, bodyFrameEndInd, ...
    wing1Res, wing1FrameStartInd, wing1FrameEndInd, ...
    wing2Res, wing2FrameStartInd, wing2FrameEndInd, mergedWingsFlag ] = ...
    hullReconstruction_mk5(params, CM_pos, all_fly_bw, body_only_bw, dlt_matrix, easyWandData,[2,1,3]);
thull = toc ;

t1 = clock ;
disp('done calculating Hulls!!!') ;
 data = hullAnalysis_mk3 (bodyRes, wing1Res, wing2Res, params, ...
    mergedWingsFlag, [], 'test');

t2 = clock ;
dt12 = t2 - t1 

% to do - measure wing vectors with respect to the vein.




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


save results flyHull  flyHullFrameStartInd  flyHullFameEndInd ...
    bodyHull   bodyHullFrameStartInd  bodyHullFameEndInd ...
    params  CM_pos  all_fly_bw body_only_bw  dlt_matrix  easyWandData ...
    cinFilenames allBG Nimages ...
    all_fly_bw_xy  body_only_bw_xy  all_fly_thresholds_xy  xcm_xy  ycm_xy  allAxlim_xy ...
    DELTA  ...
    all_fly_bw_xz  body_only_bw_xz all_fly_thresholds_xz xcm_xz ycm_xz allAxlim_xz ...
    all_fly_bw_yz  body_only_bw_yz all_fly_thresholds_yz xcm_yz ycm_yz allAxlim_yz ...
    tin tout tbody tfly 

%% segment wings on top view
[wing1_bw, wing2_bw] = findWingsOnly(all_fly_bw, body_only_bw, XY) ;