%--------------------------------------------------------------------------
% script to transfer fly data to an external hard drive. only takes
% analysis, mp4s, and README
%--------------------------------------------------------------------------
sourcePath = 'D:\Box Sync Old\VNC Sensory Lines\' ;
destinationPath = 'F:\Fly Data\Darshna Sensory Lines\' ;

% which experiments to copy over
ExprNums =  7:30 ;

% files/folders to grab
suffix_cell = {'mp4','Analysis','README.txt'} ;

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
            copyfile(sourceCurr, destCurr) ;
            disp(['Successfully copied ' sourceCurr])
        catch
            disp(['failed to copy ' sourceCurr])
            continue
        end
    end
      
end