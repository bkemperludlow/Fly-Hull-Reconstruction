% set path for experiment folder wherein analysis is still needed
pathToWatch = 'D:\Box Sync Old\VNC Motor Lines\59_25032019\' ;

% do full analysis or just angles?
justAnglesFlag = true ; 

% get experiment number from folder name--too lazy to re-enter each time
pathStruct = generatePathStruct(pathToWatch) ;
pathSplit = strsplit(pathToWatch,'\') ; 
folderSplit = strsplit(pathSplit{end-1},'_') ; 
ExprNum = str2double(folderSplit{1}) ; 

% set movie numbers that need to be analyzed
MovNum = 1:9  ; 
Nmovies = length(MovNum) ; 

% run analysis of movies
for k = 1:Nmovies
    tic
    if justAnglesFlag
        analysisPath = findMovAnalysisPath(pathStruct, MovNum(k)) ; 
        analysisPath_split = strsplit(analysisPath,'\') ; 
        savePath = strjoin({analysisPath_split{1:end-1}},'\') ; 
        estimateFlyAngles(ExprNum, MovNum(k), savePath) ; 
    else
        flyAnalysisMain(MovNum(k), ExprNum, pathStruct) ;
    end
    toc
    disp(MovNum(k)) 
end

% sort analysis folders into proper directories
moveEmptyExprFolders(pathToWatch)
