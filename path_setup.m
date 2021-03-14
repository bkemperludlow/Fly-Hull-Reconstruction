% -------------------------------------------------------------------------
% quick script to make sure proper folders are added to MATLAB for running
% Fly Hull Reconstruction code
% -------------------------------------------------------------------------
% get path to current .m file
mfile_path = matlab.desktop.editor.getActiveFilename;

% break path apart to get directory
[path, fn, ~] = fileparts(mfile_path) ; 

% get list of current folders in directory
% folderDir = dir(path) ; 
folderDir = dir(fullfile(path, '**','*')) ; 
folderDir = folderDir(3:end) ; % remove '.' and '..'
folderDir = folderDir([folderDir(:).isdir]) ; % only take folders

% define folders that we don't want to add (e.g. '.git', 'old' folders)
excludeFolders = {'.git', 'old','_old', 'under construction', ...
    'under_construction'} ; 

% check if we can exclude folders from folderDir to start
exclude_idx = arrayfun(@(x) any(contains(fullfile(x.folder,x.name), ...
    excludeFolders, 'IgnoreCase', true)) || strcmp(x.name,'.') ...
    || strcmp(x.name,'..'), folderDir) ;

folderDir = folderDir(~exclude_idx) ; 

% supress warnings for not being able to add certain paths
warning('off', 'MATLAB:mpath:privateDirectoriesNotAllowedOnPath') ;
warning('off', 'MATLAB:mpath:packageDirectoriesNotAllowedOnPath') ;

% add remaining folders
for n = 1:length(folderDir)
     % add curent path
     addpath(fullfile(folderDir(n).folder, folderDir(n).name)) ;
end

% turn warning back on
warning('on', 'MATLAB:mpath:privateDirectoriesNotAllowedOnPath') ;
warning('on', 'MATLAB:mpath:packageDirectoriesNotAllowedOnPath') ;