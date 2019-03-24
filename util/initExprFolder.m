%==========================================================================
% Script to create experiment folder in 'rootPath' directory
%==========================================================================

% this is the root folder that we want to have our expt in
rootPath = 'D:\Box Sync Old\VNC Motor Lines\' ; 
readmePath = 'D:\' ;
readmeFilename_orig = 'README-template.txt' ; 
readmeFilename_new = 'README.txt' ; 

% get index for expt folders in root
exprDir = dir(rootPath) ; 
exprDir = exprDir(3:end) ; 
exprNums = arrayfun(@(x) str2double(x.name(1:2)),exprDir) ; 

% define new folder
exprNumNew = max(exprNums) + 1 ; 
exprNumNewStr = num2str(exprNumNew, '%02d') ;

t_now = datetime ; 
t_nowStr = datestr(t_now,'ddmmyyyy') ;

newExprFolderName = strcat(exprNumNewStr,'_',t_nowStr) ; 

% make directory
newExprPath = fullfile(rootPath, newExprFolderName) ;
mkdir(newExprPath) ;
mkdir(fullfile(newExprPath,'calibration')) ;

% move over README file 
copyfile(fullfile(readmePath,readmeFilename_orig),newExprPath) ;
movefile(fullfile(newExprPath,readmeFilename_orig),...
    fullfile(newExprPath, readmeFilename_new)) ;