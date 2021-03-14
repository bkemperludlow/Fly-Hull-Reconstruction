%--------------------------------------------------------------------------
% function to transfer fly data to an external hard drive. only takes
% analysis, mp4s, and README
%
% INPUTS:
%   - sourcePath: root directory in which experiment folders are stored on
%       this machine (e.g. 'D:\Box Sync Old\Opto Silencing\')
%   - destinationPath: root directory to which files should be copied,
%       likely on external drive (e.g. 'F:\Fly Data\Opto Silencing\')
%   - ExprNums: array of experiment numbers to copy from source to
%       destination (e.g. [46, 47, 48])
%--------------------------------------------------------------------------
function copyFlyDataFun(sourcePath, destinationPath, ExprNums)
% --------------------------------
%% inputs and params
% source root directory
if ~exist('sourcePath','var') || isempty(sourcePath) 
   fprintf('Please select SOURCE directory \n')
   sourcePath = uigetdir(pwd, 'Select SOURCE directory') ; 
   fprintf('SOURCE PATH = %s \n', sourcePath)
end
% destination root directory
if ~exist('destinationPath','var') || isempty(destinationPath) 
   fprintf('Please select DESTINATION directory \n')
   destinationPath = uigetdir(pwd, 'Select DESTINATION directory') ; 
   fprintf('DESTINATION PATH = %s \n', destinationPath)
end
% experiment folder numbers
if ~exist('ExprNums','var') || isempty(ExprNums) 
   ExprNums = input(['Please enter array of experiment numbers to copy ',...
       'e.g. [45, 46, 47]']) ; 
   fprintf('Selected experiments: \n %d \n', ExprNums)
end

% ----------------------------------------------------------
% check that we got valid inputs
if ~exist(sourcePath, 'dir') 
    fprintf('Invalid source path selection -- quitting \n')
    return
end
if isempty(ExprNums) 
    fprintf('Invalid experiment number selection(s) -- quitting \n')
    return
end
%  try make destination path if it does not already exist
if ~exist(destinationPath, 'dir') 
    fprintf('Destination path does not exist -- creating now \n')
    status = mkdir(destinationPath) ;
    if ~status
        fprintf('Could not create destination path \n')
        return
    end
end
% -----------------------------------
% sort experiment numbers
ExprNums = sort(ExprNums)  ;

% files/folders to grab
suffix_cell = {'mp4','Analysis','README.txt', ...
    'calibration\calibration_dltCoefs.csv', 'calibration\wandPoints.csv'} ;

% ----------------------------------------------------------------
%% get directory structure for source and check folder existence
sourceDir = dir(sourcePath) ;
sourceDir = sourceDir(3:end) ;
sourceDir = sourceDir([sourceDir(:).isdir]) ;

% find experiment numbers
exprFolderExp = '(?<exprNum>\d+)_(?<datenum>\d+)' ;
sourceDirRegExp = arrayfun(@(x) regexp(x.name, exprFolderExp, 'names'),...
    sourceDir, 'UniformOutput', 0) ;
% remove folders from dir struct that don't conform to our experiment
% folder format
empty_idx = cellfun(@(y) isempty(y), sourceDirRegExp) ; 
sourceDirRegExp = sourceDirRegExp(~empty_idx) ; 

% get experiment numbers as double values and make sure the full list
% contains all the experiment folders we want to copy over
sourceDirExprNums = cellfun(@(y) str2double(y.exprNum), sourceDirRegExp) ; 
match_idx = ismember(ExprNums, sourceDirExprNums) ; 

if sum(match_idx) ~= length(ExprNums)
    fprintf('Warning: could not locate the following experiments: \n %d \n',...
        ExprNums(~match_idx))
end

% remove any non-matching experiment numbers
ExprNums = ExprNums(match_idx) ; 

% ---------------------------------------------------
%% loop over experiment folders
for i = 1:length(ExprNums)
    ExprNumCurr = ExprNums(i) ;
    exprInd = (sourceDirExprNums == ExprNumCurr) ;
    
    if sum(exprInd) ~= 1
        disp('Could not find experiment folder')
        continue
    end
    
    % make directory in destination
    destFolder = fullfile(destinationPath, sourceDir(exprInd).name) ;
    mkdir(destFolder) ;
    
    % copy over folders/files
    for j = 1:length(suffix_cell)
        try
            sourceCurr = fullfile(sourceDir(exprInd).folder, ...
                sourceDir(exprInd).name, suffix_cell{j}) ;
            destCurr = fullfile(destFolder, suffix_cell{j}) ;
            [status, msg] = copyfile(sourceCurr, destCurr) ;
            disp(['Successfully copied ' sourceCurr])
        catch
            disp(['failed to copy ' sourceCurr])
            continue
        end
    end
      
end
end