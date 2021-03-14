%--------------------------------------------------------------------------
% generates a structure containing relevant file paths for a given
% experiment. 
%
% pathToWatch is the root folder for the experiment, e.g. :
%       'D:\Box Sync Old\VNC Motor Lines\11_25032018\'
%--------------------------------------------------------------------------
function pathStruct = generatePathStruct(pathToWatch)
% -------------------------------------------------------------------------
%% add the various paths to the structure
pathStruct = struct ; 
pathStruct.root = pathToWatch ; 
pathStruct.save = fullfile(pathToWatch, 'Analysis') ; 
pathStruct.calibration = fullfile(pathToWatch, 'calibration') ; 
pathStruct.errorLog = fullfile(pathToWatch, 'errorlog.txt') ;
pathStruct.possibleFT = fullfile(pathToWatch, 'Possible False Triggers\') ; 
pathStruct.pitchUp = fullfile(pathToWatch, 'Analysis\Pitch Up\') ; 
pathStruct.pitchDown = fullfile(pathToWatch, 'Analysis\Pitch Down\') ; 
pathStruct.rollRight = fullfile(pathToWatch, 'Analysis\Roll Right\') ; 
pathStruct.rollLeft = fullfile(pathToWatch, 'Analysis\Roll Left\') ; 
pathStruct.noPert = fullfile(pathToWatch, 'Analysis\No Perturbation\') ; 
pathStruct.probNoPert = fullfile(pathToWatch, 'Analysis\Probably No Perturbation\') ; 
pathStruct.undeterminedPert = fullfile(pathToWatch, 'Analysis\Undetermined Perturbation\') ;
pathStruct.other = fullfile(pathToWatch, 'Analysis\Other\') ; 
pathStruct.unsorted = fullfile(pathToWatch, 'Analysis\Unsorted\') ; 
pathStruct.mp4 = fullfile(pathToWatch, 'mp4') ;

% -------------------------------------------------------------------------
%% get experiment number 
pathToWatch_split = strsplit(pathToWatch, '\') ;
emptyInd = cellfun(@(y) isempty(y), pathToWatch_split) ; 
pathToWatch_split = pathToWatch_split(~emptyInd) ; 

exprFolder = pathToWatch_split{end} ; 
exprFolder_split = strsplit(exprFolder, '_') ; 
exprNumStr = exprFolder_split{1} ; 

ExprNum = str2double(exprNumStr) ; 
pathStruct.ExprNum = ExprNum ; 

% -------------------------------------------------------------------------
%% create the directory folders if they don't already exist
pathStructFields = {'root', 'save', 'possibleFT', 'pitchUp', 'pitchDown',...
    'rollRight', 'rollLeft', 'noPert', 'probNoPert','undeterminedPert',...
    'unsorted','other','mp4'} ; 
for i = 1:length(pathStructFields)
    if ~exist(pathStruct.(pathStructFields{i}), 'dir')
        mkdir(pathStruct.(pathStructFields{i}))
    end
end

end