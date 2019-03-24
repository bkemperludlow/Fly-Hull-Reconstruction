function [ExprNum, MovNum, tinArray, toutArray, pertFlags, TwoFliesXY, TwoFliesXZ, TwoFliesYZ] = ...
    generateAnalysisList(movieSummaryPath)
%This should output the arrays ExprNum, MovNum, tinArray, toutArray,
%pertFlags, TwoFliesXY, TwoFliesXZ, and TwoFliesYZ for the runManyAnalyses
%script

%movieSummaryPath = 'G:\Janelia Flies\kir2.1 flies round 2\' ; 
movieSummary = importdata([movieSummaryPath 'Movie Summary.xlsx']) ;

movieSummaryData = movieSummary.data ; 
movieSummaryTextData = movieSummary.textdata ;

ExprNum = [] ; 
MovNum = [] ; 
tinArray = [] ; 
toutArray = [] ; 
pertFlags = [] ; 

cc = 1 ; 

for i = [1 2 5 6 7] %i just know that there are 5 experiments I want to look at; should adapt this later
    
    sheetNameStr = ['Exp' num2str(i)] ; 
    dataCurr = movieSummaryData.(sheetNameStr) ;
    textDataCurr = movieSummaryTextData.(sheetNameStr) ;
    
    Nmovies = size(dataCurr,1) ; 
    for j = 1:Nmovies
        
        pertTypeCurr = textDataCurr{2+j, 2} ; 
        pertDirectionCurr = textDataCurr{2+j, 3} ; 
        
        if strcmp(pertTypeCurr,'FT')
            continue ; 
        end
        
        ExprNum(cc) = i ; 
        MovNum(cc) = dataCurr(j,1) ; 
        tinArray(cc) = 8*dataCurr(j, 5) ; %convert from millisecond to frame number
        toutArray(cc) = 8*dataCurr(j, 6) ; 
        
        if strcmp(pertTypeCurr,'Pitch') && strcmp(pertDirectionCurr,'Up')
            pertFlags(cc) = 1 ; 
        elseif strcmp(pertTypeCurr,'Pitch') && strcmp(pertDirectionCurr,'Down')
            pertFlags(cc) = -1 ; 
        elseif strcmp(pertTypeCurr,'Roll') && strcmp(pertDirectionCurr,'Right')
            pertFlags(cc) = 2 ; 
        elseif strcmp(pertTypeCurr,'Roll') && strcmp(pertDirectionCurr,'Left')
            pertFlags(cc) = -2 ; 
        elseif strcmp(pertTypeCurr,'None') 
            pertFlags(cc) = 0 ; 
        else 
            pertFlags(cc) = 3 ; 
        end
        
        cc = cc + 1 ; 
       
    end   

end

TwoFliesXY = zeros(1, cc - 1) ; 
TwoFliesXZ = zeros(1, cc - 1) ; 
TwoFliesYZ = zeros(1, cc - 1) ; 

