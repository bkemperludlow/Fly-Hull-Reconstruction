%--------------------------------------------------------------------------
% function used in main analysis to estimate roll vectors and calculate
% body/wing euler angles
%--------------------------------------------------------------------------
function data = estimateFlyAngles(exprNum, movNum, savePath,...
    largePertFlag, overWriteFlag)
%--------------------------------------------------------------------------
%% params and inputs
if ~exist('largePertFlag','var') || isempty(largePertFlag)
    largePertFlag = false ;
end
if ~exist('overWriteFlag','var')
    overWriteFlag = true ;
end

%--------------------------------------------------------------------------
%% get file and directory info
movNumStr = num2str(movNum,'%03d') ; 
folderName = ['Expr_' num2str(exprNum) '_mov_' movNumStr] ;
datapath = fullfile(savePath, folderName) ;

datafilename_test = fullfile(datapath, [folderName '_test.mat']) ;
%{
datafilename_cleaned = fullfile(datapath, [folderName '_cleaned.mat']) ;
datafilename_refined = fullfile(datapath, [folderName '_iterRefine.mat']) ;
datafilename_corrected = fullfile(datapath,['Expr' num2str(exprNum) 'mov' movNumStr,'_Data_manually_corrected.mat']);
if exist(datafilename_corrected,'file')
    datafilename = datafilename_corrected;
elseif exist(datafilename_refined, 'file')
    datafilename = datafilename_refined;
elseif exist(datafilename_cleaned,'file')
    datafilename = datafilename_cleaned ;
else 
%}
    datafilename = datafilename_test ;
%end
    
%cd(datapath)

if (exist(datafilename,'file') == 2) && ~overWriteFlag
    % if we're not supposed to overwrite
   data = importdata(datafilename) ; 
   return
elseif (exist(datafilename,'file') == 2) && overWriteFlag
    % calculate angles
    data = calcAnglesMain(datafilename, largePertFlag, true, true) ;
   
   % close plot windows
   close all
else 
   analysis_results = load(fullfile(savePath, folderName, ...
       [folderName '_results.mat'])) ; 
   set(0, 'ShowHiddenHandles', 'on')
   try
       close 'easyWand 5'
   catch
       disp('No easy wand window')
   end
    
   % estimate roll vectors
   data = analysis_results.data ; 
   [rhoTimes, rollVectors] = estimateRollVector(data,largePertFlag) ;
   data.rhoTimes = rhoTimes ;
   data.rollVectors = rollVectors ;
   save(datafilename, 'data')
   
   % calculate angles
   data = calcAnglesMain(datafilename, largePertFlag, true, true) ;
   
   % close plot windows
   close all
end
end
%end