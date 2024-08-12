% -------------------------------------------------------------------------
% to be used alongside "checkFalseTriggers.m" ; runs "moveToFalseTrig" for
% a set of movies using hacky workaround
% -------------------------------------------------------------------------
rootPath = 'D:\Box Sync Old\VNC MN Chrimson\84_24102021\' ;
falseTrigStruct = importdata(fullfile(rootPath, 'falseTrigStruct.mat')) ; 
pathStruct = generatePathStruct(rootPath) ; 

for k = 1:length(falseTrigStruct)
    moveToFalseTrig(falseTrigStruct(k).MovNum, pathStruct)
    fprintf('Completed moving %d/%d movies \n', k, length(falseTrigStruct))
end