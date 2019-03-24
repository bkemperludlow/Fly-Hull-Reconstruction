rootPath = 'G:\Janelia Flies\kir2.1 flies\Analysis\' ;
%cd(rootPath)

pitchUp_ExprNum = [3 3 3 4 4 4 4 4 4 8 8 8 8 10 10 10 11 11 12 12 13 13 13 13] ; 
pitchUp_MovNum = [10 30 37 4 10 106 137 147 202 13 23 29 96 38 43 45 20 38 11 64 103 130 134 138] ;
pitchDown_ExprNum = [3 3 3 3 4 4 6 10 7 7 7 8 8 10 10 11 12 12] ;
pitchDown_MovNum = [25 48 96 124 68 109 98 19 53 102 145 16 70 47 51 9 23 57] ;

halfStrokeStruct = struct('ExprNum', [], 'MovNum', [], 'fwdStrokes',[],'backStrokes',[],...
    'pitchType',[],'flyType',[]);

