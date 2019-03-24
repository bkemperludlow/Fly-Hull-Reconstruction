%function [xy_cm, xz_cm, yz_cm] = findImageCM(ExprNum, MovNum,data) 

%% Gives center of mass in normal image coordinates. 
%Need to convert to real coordinates for animation

%--------------------------------------------------------------------------
%% Set parameters

ExprNum = 7;
MovNum = 8 ;

%% Load data
phantomSDK_setPath ;
LoadPhantomLibraries ;

XZ = 2 ;
XY = 3 ;
YZ = 1 ;

if MovNum < 10
    zstr = '00' ;
elseif MovNum < 100
    zstr = '0' ;
else
    zstr = '' ;
end

dataPath = 'F:\luca\07_230514\cines' ;
cinFilenames = cell(1,3) ;
cinFilenames{XY} = [dataPath strcat('\xy_',zstr,num2str(MovNum),'.cin')] ;
cinFilenames{XZ} = [dataPath strcat('\xz_',zstr,num2str(MovNum),'.cin')] ;
cinFilenames{YZ} = [dataPath strcat('\yz_',zstr,num2str(MovNum),'.cin')] ;
    
tin = data.params.startTrackingTime ;
tout = data.params.endTrackingTime ;

%% Find backgrounds of images (taken from 'beginning.m')
% background images in two formats (matrix and cell). I probably only use
% the matrix form.
allBG   = zeros(3, 512, 512,'uint8') ; % cell(1,3) ;
allBGcell = cell(1,3) ;

% estimations for the time the fly comes in and out of the FOV of each
% camera. this part of the automation can be improved. check it or just set
% the "tin" and "tout" manually later.
allTin  = zeros(3,1) ;
allTout = zeros(3,1) ;

tic

disp('Finding background images. Batch mode...') ;

job1 = batch('findBG',3,{cinFilenames{1}}) ;
job2 = batch('findBG',3,{cinFilenames{2}}) ;
job3 = batch('findBG',3,{cinFilenames{3}}) ;

wait(job1) ;
wait(job2) ;
wait(job3) ;

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

clear out1 out2 out3
batchtime = toc ;

disp('Done finding background images...') ;
disp(' ') ;
%{
tic ;
for cam = 1:3 
    [bg, tin_curr, tout_curr] = findBG(cinFilenames{cam}); % finds background here.
    allBG(cam,:,:) = bg ;
    allBGcell{cam} = bg ;
    allTin(cam) = tin_curr ;
    allTout(cam) = tout_curr ;    
end
regtime = toc 
%}

%% Find centers of mass using binaryThreshold.m

tic 

cam = XY ;
%jobXY = batch('binaryThreshold',7,{squeeze(allBG(cam,:,:)), cinFilenames{cam}, tin, tout}) ;
[~, ~, ~, xcm_xy, ycm_xy, ~, ~] = ...
    binaryThreshold( squeeze(allBG(cam,:,:)) , cinFilenames{cam}, tin, tout) ;

cam = XZ ;
%jobXZ = batch('binaryThreshold',7,{squeeze(allBG(cam,:,:)), cinFilenames{cam}, tin, tout}) ;
[~, ~, ~, xcm_xz, ycm_xz, ~, ~] = ...
    binaryThreshold( squeeze(allBG(cam,:,:)) , cinFilenames{cam}, tin, tout) ;


cam = YZ ;
%jobYZ = batch('binaryThreshold',7,{squeeze(allBG(cam,:,:)), cinFilenames{cam}, tin, tout}) ;
[~, ~, ~, xcm_yz, ycm_yz, ~, ~] = ...
    binaryThreshold( squeeze(allBG(cam,:,:)) , cinFilenames{cam}, tin, tout) ;

%% Store the output

xy_cm = [xcm_xy, ycm_xy] ;
xz_cm = [xcm_xz, ycm_xz] ;
yz_cm = [xcm_yz, ycm_yz] ; 
