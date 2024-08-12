%--------------------------------------------------------------------------
% function to load in data for a controller fit. kind of old, could use
% some improvements
%--------------------------------------------------------------------------

function [data, dataPath, pertType, ExprNum, dataFilename] = ...
    loadPertDataStruct(rootPath, MovNum)
%% get folder and potential data file names 
pathStruct = generatePathStruct(rootPath) ; 
dataPath = findMovAnalysisPath(pathStruct, MovNum) ; 
ExprNum = pathStruct.ExprNum ; 

MovNumStr = num2str(MovNum,'%03d') ; 
dataFilename1 = ['Expr' num2str(ExprNum) 'mov' MovNumStr ...
    '_Data_manually_corrected.mat'] ;
dataFilename2 = ['Expr' num2str(ExprNum) 'mov' MovNumStr ...
    '_Data_manually_corrected_onlyPhi.mat'] ;
dataFilename3 = ['Expr' num2str(ExprNum) 'mov' MovNumStr ...
    '_Data_manually_corrected_onlyPhiFront.mat'] ;
dataFilename4 = ['Expr_' num2str(ExprNum) '_mov_' MovNumStr '_cleaned.mat'] ;
dataFilename5 = ['Expr_' num2str(ExprNum) '_mov_' MovNumStr '_test.mat'] ;

%% get data folder
if contains(dataPath,'Pitch Up','IgnoreCase',true)
    pertType = 1 ; 
elseif contains(dataPath,'Pitch Down','IgnoreCase',true)
    pertType = -1 ; 
elseif contains(dataPath,'Roll Right','IgnoreCase',true)
    pertType = 2 ; 
elseif contains(dataPath,'Roll Left','IgnoreCase',true)
    pertType = -2 ;
elseif contains(dataPath,'Yaw Right','IgnoreCase',true)
    pertType = 3 ;
elseif contains(dataPath,'Yaw Left','IgnoreCase',true)
    pertType = -3 ;
else
    pertType = 0 ; 
end

%% try to load data. first try to find a manually corrected version
if (exist(fullfile(dataPath,dataFilename1),'file') == 2)
    data = importdata(fullfile(dataPath,dataFilename1)) ;
    dataFilename = dataFilename1 ; 
    
elseif (exist(fullfile(dataPath,dataFilename2),'file') == 2)
    data = importdata(fullfile(dataPath,dataFilename2)) ;
    dataFilename = dataFilename2 ; 
    
elseif (exist(fullfile(dataPath,dataFilename3),'file') == 2)
    data = importdata(fullfile(dataPath,dataFilename3)) ;
    dataFilename = dataFilename3 ; 
    
elseif (exist(fullfile(dataPath,dataFilename4),'file') == 2)
    data = importdata(fullfile(dataPath,dataFilename4)) ;
    dataFilename = dataFilename4 ; 
    
elseif (exist(fullfile(dataPath,dataFilename5),'file') == 2)
    data = importdata(fullfile(dataPath,dataFilename5)) ;
    dataFilename = dataFilename5 ; 
    
else
    disp('No data file')
    data = [] ;
    return
end



end