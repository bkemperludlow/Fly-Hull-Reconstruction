% -------------------------------------------------------------------------
% function to compress all cine files in a given experiment folder
% using Cine2SparseSam.m (by experiment folder, we mean folder of the form
% XX_DDMMYYY, as is typical of flight experiments)
% 
% ** RUN AS "CineToSparseFormat.compressAllMoviesInFolder(argin)" since it
% is part of a class
%  
% NB: will maybe eventually want to add something to move raw cine files;
% however, that could probably just go in another script
% -------------------------------------------------------------------------
function [] = compressAllMoviesInFolder(exprRoot, bgMethod, ...
    defaultFrameRate, parallelFlag)
% ---------------------
%% 
% inputs and params
if ~exist('exprRoot','var') || isempty(exprRoot)
    startPath = 'D:\Box Sync Old\VNC MN Chrimson\81_14082021\' ; % EDIT FOR CONVENIENCE
    exprRoot = uigetdir(startPath, ['Select experiment folder for' ...
        'cine compression']) ; 
    if isequal(exprRoot,0)
        disp('User selected Cancel; end of session')
        return
    end
end
if ~exist('bgMethod','var') || isempty(bgMethod)
    bgMethod='Single movie max better';
     %{
    'Single movie max better'
    'Single movie max'
    'Single movie mean'
    'Single movie median'
    'Manual input'
    %}
end
if ~exist('defaultFrameRate','var') || isempty(defaultFrameRate)
    % default frame rate for cine (in case we can't get it from metadata)
    defaultFrameRate = 8000 ; % frames per second
end
if ~exist('parallelFlag','var') || isempty(parallelFlag)
    % run compression of multiple cines in parallel?
    parallelFlag = false ;
end

% --------------
% misc. params
validExts = {'.cin', '.cine'} ; % phantom camera video file types
checkBGFlag = false ; % can't do this in real time for batch

% -------------------------
%% 
% get all phantom camera video files in exprRoot
exprDir = cellfun(@(y) dir(fullfile(exprRoot, ['*' y])), validExts, ...
    'UniformOutput',0) ; 
exprDir = cat(1, exprDir{:}) ; 

% get cin/cine filenames from directory structure
filenames = arrayfun(@(x) fullfile(x.folder, x.name), exprDir, ...
    'UniformOutput',0) ; 

% sort filenames by save time, so (ideally) movies will be nicely grouped
% in file explorer
datenums = [exprDir.datenum] ; 
[~, sort_ind] = sort(datenums) ; 
filenames = filenames(sort_ind) ; 
% -------------------------
%%
% run compression on all video files 
%  (loop happens inside function, unless we want to parallelize)
if parallelFlag
    % execute parallel for loop
    parfor k = 1:length(filenames)
        % run compression on individual files
        CineToSparseFormat.Cine2SparseSam(filenames{k}, bgMethod, ...
            checkBGFlag, defaultFrameRate) ;
    end
    
    % shut down parallel pool
    delete(gcp) ; 

else
    % otherwise just run script as normal
    CineToSparseFormat.Cine2SparseSam(filenames, bgMethod, ...
        checkBGFlag, defaultFrameRate) ;
end

% print to signify completion
fprintf('Completed compressing movie files in %s \n', exprRoot)

end