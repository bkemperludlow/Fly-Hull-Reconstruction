% -------------------------------------------------------------------------
% function to redo hullAnalysis (mostly to account for errors in long body
% axis calculation)
% -------------------------------------------------------------------------
function data = redoHullAnalysis(analysisOutput)
% params
plotHullFlag = false ;
saveHullFigFlag = false ;
hullFigPath = [] ;

% read in data
bodyRes = analysisOutput.bodyRes ;
wing1Res = analysisOutput.wing1Res ;
wing2Res = analysisOutput.wing2Res ;
params = analysisOutput.params ; 
mergedWingsFlag = analysisOutput.mergedWingsFlag ; 

% perform hull analysis
data = hullAnalysis_mk5 (bodyRes, wing1Res, wing2Res, params, ...
    mergedWingsFlag, [], 'test', plotHullFlag, saveHullFigFlag, ...
    hullFigPath);

end

