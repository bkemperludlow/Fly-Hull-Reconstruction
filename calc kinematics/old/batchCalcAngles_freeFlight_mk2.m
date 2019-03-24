%mk2: Changed to function work better with the automated analysis
%Can still use original for batch?

function data = batchCalcAngles_freeFlight_mk2(exprNum, movNum, savePath, bigPitchFlag)
%rootPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Roll Right' ;
%rootPath = 'D:\Raymond Analysis Test\camera test automation\Analysis' ;

if nargin < 4
    bigPitchFlag = false ;
end
cd(savePath) ;
overWriteFlag = true ; 

%cd([savePath '\No Perturbation'])
%movieDir = dir(prefixStr) ;
movNumStr = num2str(movNum,'%03d') ; 
movieDir = dir(['Expr_' num2str(exprNum) '_mov_' movNumStr]) ;
%Nmovies = length(movieDir) ;

%for i = 1:length(movieDir)
   %folderName = movieDir.name ; 
folderName = ['Expr_' num2str(exprNum) '_mov_' movNumStr] ;
datapath = [savePath '\' folderName '\'] ;
exprNum = str2double(folderName(end-9:end-8)) ;
if isnan(exprNum)
   exprNum = str2double(folderName(end-8)) ;
end
movNum = str2double(folderName(end-2:end)) ;

datafilename = [datapath '\' folderName '_test.mat'] ;

cd(datapath)

if (exist(datafilename,'file') == 2) && ~overWriteFlag
   data = importdata(datafilename) ; 
   return
else 
   load([savePath '\' folderName '\' folderName '_results.mat'])
   set(0, 'ShowHiddenHandles', 'on')
   try
       close 'easyWand 5'
   catch
       disp('No easy wand window')
   end

   [rhoTimes, rollVectors] = estimateRollVector_v2(data) ;
   data.rhoTimes = rhoTimes ;
   data.rollVectors = rollVectors ;

   run quick_and_dirty
   %run quick_and_dirty_auto
   
   save(datafilename, 'data')
   %keyboard ;
   close all
end
end
%end