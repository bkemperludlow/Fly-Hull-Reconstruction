% -------------------------------------------------------------------------
% temporary script to run "cleanUpAnalysesFun.m" for several different
% experiments
% -------------------------------------------------------------------------
clustFlag = false ; 
largePertFlag = true ; 

% -----------------------------------------------------------------

pathToWatch = 'D:\Fly Data\VNC MN Chrimson\60_18102020\' ;
MovNum = [10, 16, 39, 42, 88, 92:134] ; 

analysisType = 'clean_wings' ; 
cleanUpAnalysesFun(pathToWatch, analysisType, MovNum, clustFlag, ...
    largePertFlag) ; 

analysisType = 'extreme_roll' ; 
cleanUpAnalysesFun(pathToWatch, analysisType, MovNum, clustFlag, ...
    largePertFlag) ; 

% -----------------------------------------------------------------

