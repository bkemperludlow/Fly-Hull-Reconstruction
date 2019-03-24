% set path for experiment folder wherein analysis is still needed
pathToWatch = 'D:\Box Sync Old\VNC Motor Lines\56_23022019\' ;

% get experiment number from folder name--too lazy to re-enter each time
pathStruct = generatePathStruct(pathToWatch) ;
pathSplit = strsplit(pathToWatch,'\') ; 
folderSplit = strsplit(pathSplit{end-1},'_') ; 
ExprNum = str2double(folderSplit{1}) ; 

% set movie numbers that need to be analyzed
MovNum = []  ; 
Nmovies = length(MovNum) ; 

% run analysis of movies
for k = 1:Nmovies
    tic
    flyAnalysisMain(MovNum(k), ExprNum, pathStruct) ;
    toc
    disp(MovNum(k)) 
end

% sort analysis folders into proper directories
moveEmptyExprFolders(pathToWatch)
