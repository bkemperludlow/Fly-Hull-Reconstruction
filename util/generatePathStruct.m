function pathStruct = generatePathStruct(pathToWatch)

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

mkdir(pathStruct.root) ; 
mkdir(pathStruct.save) ; 
mkdir(pathStruct.possibleFT) ; 
mkdir(pathStruct.pitchUp ) ; 
mkdir(pathStruct.pitchDown) ; 
mkdir(pathStruct.rollRight) ; 
mkdir(pathStruct.rollLeft) ; 
mkdir(pathStruct.noPert) ;
mkdir(pathStruct.probNoPert) ;
mkdir(pathStruct.undeterminedPert) ;
mkdir(pathStruct.unsorted) ;
mkdir(pathStruct.other) ;
mkdir(pathStruct.mp4) ;


end