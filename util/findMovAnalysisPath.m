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
FAIL_IND = arrayfun(@(x) contains(x.folder, 'Failed Analysis'), ...
    analysisDir) ; 
IND = IND & ~FAIL_IND ; 

if (sum(IND) < 1)
   fprintf('Could not find folder for movie %d \n', MovNum)
   movAnalysisPath = [] ; 
   return
elseif (sum(IND) > 1)
    fprintf('Multiple folders for movie %d -- picking one \n', MovNum)
    
    % loop through the different folder options and find memory stored in
    % each 
    candidate_idx = find(IND) ; 
    folderMem = zeros(length(candidate_idx),1) ; 
    for i = 1:length(candidate_idx)
       idx_curr = candidate_idx(i) ; 
       dir_curr = dir(fullfile(analysisDir(idx_curr).folder, ...
           analysisDir(idx_curr).name)) ; 
       folderMem(i) = sum([dir_curr.bytes]) ; 
    end
    
    % take the folder with the most occupied memory
    [~, max_ind] = max(folderMem) ;
    IND_NEW = false(size(IND)) ; 
    IND_NEW(candidate_idx(max_ind)) = true ; 
    
    % upade index to new folder 
    IND = IND_NEW ; 
end

% get path
movAnalysisPath = fullfile(analysisDir(IND).folder, analysisDir(IND).name) ; 

end