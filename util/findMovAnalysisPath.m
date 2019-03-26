%--------------------------------------------------------------------------
% function to find the path to a movie's analysis folder
%
% for example, if we have a data file whose full path is:
%   'D:\Box Sync Old\VNC MN Chrimson\25_08122018\Analysis\Unsorted\Expr_25_mov_000\Expr_25_mov_000_test.mat'
% this function should return:
%   'D:\Box Sync Old\VNC MN Chrimson\25_08122018\Analysis\Unsorted\Expr_25_mov_000'
%--------------------------------------------------------------------------
function movAnalysisPath = findMovAnalysisPath(pathStruct, MovNum)
% get directory structure within umbrella Analysis folder
analysisDir = dir([pathStruct.save '\**\Expr*']) ;
analysisDir = analysisDir([analysisDir(:).isdir]) ;

% make sure that the folders we're looking at correspond to individual
% analysis folders (i.e. have the format we expect)
folderSplit = arrayfun(@(x) strsplit(x.name,'_'),analysisDir, ...
    'UniformOutput', 0) ; 
movFolderInd = cellfun(@(y) (length(y)==4), folderSplit) ; 

analysisDir = analysisDir(movFolderInd) ; 
folderSplit = folderSplit(movFolderInd) ; 

% find folder with correct movie number
MovNumStr = num2str(MovNum, '%03d') ; 
IND = cellfun(@(y) strcmp(y{4}, MovNumStr), folderSplit) ;

if (sum(IND) ~= 1)
   disp('could not find unique folder for this movie')
   movAnalysisPath = [] ; 
   return
end

% get path
movAnalysisPath = fullfile(analysisDir(IND).folder, analysisDir(IND).name) ; 

end