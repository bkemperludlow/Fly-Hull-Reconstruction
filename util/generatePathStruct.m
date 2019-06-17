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
pathStruct.save = [pathToWatch 'Analysis'] ; 
pathStruct.calibration = [pathToWatch 'calibration'] ; 
pathStruct.errorLog = [pathToWatch 'errorlog.txt'] ;
pathStruct.possibleFT = [pathToWatch 'Possible False Triggers\'] ; 
pathStruct.pitchUp = [pathToWatch 'Analysis\Pitch Up\'] ; 
pathStruct.pitchDown = [pathToWatch 'Analysis\Pitch Down\'] ; 
pathStruct.rollRight = [pathToWatch 'Analysis\Roll Right\'] ; 
pathStruct.rollLeft = [pathToWatch 'Analysis\Roll Left\'] ; 
pathStruct.noPert = [pathToWatch 'Analysis\No Perturbation\'] ; 
pathStruct.probNoPert = [pathToWatch 'Analysis\Probably No Perturbation\'] ; 
pathStruct.undeterminedPert = [pathToWatch 'Analysis\Undetermined Perturbation\'] ;
pathStruct.other = [pathToWatch 'Analysis\Other\'] ; 
pathStruct.unsorted = [pathToWatch 'Analysis\Unsorted\'] ; 
pathStruct.mp4 = [pathToWatch 'mp4'] ;

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