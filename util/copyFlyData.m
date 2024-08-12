%--------------------------------------------------------------------------
% script to transfer fly data to an external hard drive. only takes
% analysis, mp4s, and README
%--------------------------------------------------------------------------

sourcePath = 'D:\VNC Motor Lines\' ; %'D:\Box Sync Old\VNC Motor Lines\' ;
destinationPath = 'C:\Users\Kemper\Documents\Fly Experiment Data\VNC Motor Lines\' ; % 'F:\Fly Data\VNC Motor Lines\' ;

% which experiments to copy over
ExprNums = [24:30,33:37,64,65,71:73,77,79:86,102:111,126:129,133:151]  ;


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
    
%     % also make calibration directory, if using
%     if any(cellfun(@(y) contains(y,'calibration'), suffix_cell))
%        mkdir(fullfile(destFolder,'calibration')) 
%     end
    
    % copy over folders/files
    for j = 1:length(suffix_cell)
        % make parent directories, if needed
        suffix_split = strsplit(suffix_cell{j},'\') ;
        if length(suffix_split) > 1
            newFolder = fullfile(destFolder, ...
                strjoin(suffix_split(1:end-1),'\')) ;
            if ~exist(newFolder,'dir')
               mkdir(newFolder) 
            end
        end
        
        % try to transfer data
        sourceCurr = fullfile(sourceDir(exprInd).folder, ...
            sourceDir(exprInd).name, suffix_cell{j}) ;
        destCurr = fullfile(destFolder, suffix_cell{j}) ;
        [status, msg] = copyfile(sourceCurr, destCurr) ;
        
        % print out whether or not it worked
        if status
            fprintf('Successfully copied %s \n', sourceCurr)
        else
            fprintf('Failed to copy %s \n', sourceCurr)
            disp(msg)
        end
    end    
end