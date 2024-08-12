% -------------------------------------------------
% script to perform various analysis tasks
% -------------------------------------------------
% clean up analysis of Opto Silencing / Expr 89
pathToWatch = 'D:\Box Sync Old\Opto Silencing\89_29102021\' ; 
analysisType = 'clean_wings' ; 
MovNum = 1:26 ; 
clustFlag = false ;
largePertFlag = false; 

cleanUpAnalysesFun(pathToWatch, analysisType, MovNum, clustFlag, ...
    largePertFlag) ;



