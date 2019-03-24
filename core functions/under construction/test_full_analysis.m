%--------------------------------------------------------------------------
% quick and dirty script to make sure i haven't fucked up the code by
% modifying it
%--------------------------------------------------------------------------

%% first normal
rootPath = 'G:\Janelia Flies\kir2.1 flies round 2\15_08052016\' ; 
clustFlag = false ;

pathStruct = generatePathStruct(rootPath) ; 
ExprNum = 15 ; 
MovNum = 4 ; 

data = flyAnalysisMain(MovNum, ExprNum, pathStruct, clustFlag) ;

%% then using the cluster 