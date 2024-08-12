%% --------------------------------------------------------------------
% function to move a set of files to "Possible False Triggers" folder
function [] = moveToFalseTrig(movNum, pathStruct)
% initialize folder in "Possible False Triggers"
newDirName = fullfile(pathStruct.possibleFT, ...
    sprintf('Expr_%d_mov_%03d',pathStruct.ExprNum, movNum)) ;
mkdir(newDirName) ;

% get files to move
cinDir =  dir(fullfile(pathStruct.root, sprintf('*_%03d.cin*',movNum))) ;
xmlDir =  dir(fullfile(pathStruct.root, sprintf('*_%03d.xml',movNum))) ;
combDir = [cinDir ; xmlDir] ;

% loop over files and move from cine root to "Possible False Triggers"
for k = 1:length(combDir)
    % current source and destination paths
    sourceCurr = fullfile(combDir(k).folder, combDir(k).name) ;
    destCurr = fullfile(newDirName, combDir(k).name) ;
   
    % try to move file -- if it fails, pause (for now)
    [status, msg] = movefile(sourceCurr, destCurr, 'f') ;
    if ~status
        disp(cmdout)
        keyboard
    end
end
end