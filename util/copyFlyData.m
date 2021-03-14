%--------------------------------------------------------------------------
% script to transfer fly data to an external hard drive. only takes
% analysis, mp4s, and README
%--------------------------------------------------------------------------
sourcePath = 'D:\Box Sync Old\Opto Silencing\' ; %'D:\Box Sync Old\VNC Motor Lines\' ;
destinationPath = 'F:\Fly Data\Opto Silencing\' ; % 'F:\Fly Data\VNC Motor Lines\' ;

% which experiments to copy over
ExprNums = sort([45])  ;

% files/folders to grab
suffix_cell = {'mp4','Analysis','README.txt', ...
    'calibration\calibration_dltCoefs.csv', 'calibration\wandPoints.csv'} ;

% get directory structure for source
sourceDir = dir(sourcePath) ;
sourceDir = sourceDir(3:end) ;
sourceDir = sourceDir([sourceDir(:).isdir]) ;
% find experiment numbers
sourceDirExprNums = arrayfun(@(x) str2double(x.name(1:2)), sourceDir) ;

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