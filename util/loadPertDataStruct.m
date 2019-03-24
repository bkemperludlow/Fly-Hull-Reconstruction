%--------------------------------------------------------------------------
% function to load in data for a controller fit. kind of old, could use
% some improvements
%--------------------------------------------------------------------------

function [data, dataPath] = loadPertDataStruct(rootPath, ExprNum, MovNum, pertType)
%% get folder and potential data file names 
MovNumStr = num2str(MovNum,'%03d') ; 
folderName = ['Expr_' num2str(ExprNum) '_mov_' MovNumStr ] ;
dataFilename1 = ['Expr' num2str(ExprNum) 'mov' MovNumStr ...
    '_Data_manually_corrected.mat'] ;
dataFilename2 = ['Expr' num2str(ExprNum) 'mov' MovNumStr ...
    '_Data_manually_corrected_onlyPhi.mat'] ;
dataFilename3 = ['Expr' num2str(ExprNum) 'mov' MovNumStr ...
    '_Data_manually_corrected_onlyPhiFront.mat'] ;
dataFilename4 = ['Expr_' num2str(ExprNum) '_mov_' MovNumStr '_test.mat'] ;

%% get data folder
if pertType == 1
    dataPath = fullfile(rootPath, 'Pitch Up', folderName) ;
elseif pertType == -1 
    dataPath = fullfile(rootPath, 'Pitch Down', folderName) ;
elseif pertType == 2 
    dataPath = fullfile(rootPath, 'Roll Right', folderName) ;
elseif pertType == -2 
    dataPath = fullfile(rootPath, 'Roll Left', folderName) ; 
elseif isnan(pertType)
    dataPath = fullfile(rootPath, 'Unsorted', folderName) ;
else
    disp('Pert type not correct. Please check')
    keyboard ;
end

%% try to load data. first try to find a manually corrected version
if exist(fullfile(dataPath,dataFilename1),'file') == 2
    data = importdata(fullfile(dataPath,dataFilename1)) ;
elseif exist(fullfile(dataPath,dataFilename2),'file') == 2
    data = importdata(fullfile(dataPath,dataFilename2)) ;
elseif exist(fullfile(dataPath,dataFilename3),'file') == 2
    data = importdata(fullfile(dataPath,dataFilename3)) ;
elseif exist(fullfile(dataPath,dataFilename4),'file') == 2
    data = importdata(fullfile(dataPath,dataFilename4)) ;
else
    disp('No data file')
    keyboard ;
end



end