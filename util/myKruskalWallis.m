cd('G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Controller Analysis')
load controllerAnalysisStruct_LM

%find indices for relevant things to plot
MB204B_w1118_Ind = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    ([controllerAnalysisStruct(:).ExprNum] == 2 | ...
    [controllerAnalysisStruct(:).ExprNum] == 7)) ;

MB204B_kir_Ind = find([controllerAnalysisStruct(:).flyType] == 1 & ...
    ([controllerAnalysisStruct(:).ExprNum] == 1 | ...
    [controllerAnalysisStruct(:).ExprNum] == 8)) ;

EmptyVector_Ind = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    ([controllerAnalysisStruct(:).ExprNum] == 6 | ...
    [controllerAnalysisStruct(:).ExprNum] == 14 | [controllerAnalysisStruct(:).ExprNum] == 15 )) ;

w1118_Ind = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).ExprNum] == 5) ;

flyLine = cell(length(controllerAnalysisStruct),1) ; 
flyLine(MB204B_w1118_Ind) = {'MB204B > w1118'} ;
flyLine(MB204B_kir_Ind) = {'MB204B > kir2.1'} ;
flyLine(w1118_Ind) = {'w1118 > kir2.1'} ;
flyLine(EmptyVector_Ind) = {'EmptyVector > kir2.1'} ;

K_i = [controllerAnalysisStruct(:).K_i] ; 
K_p = [controllerAnalysisStruct(:).K_p] ; 
gainRatio = K_i ./ K_p ;
deltaT = [controllerAnalysisStruct(:).deltaT] ; 
rms = [controllerAnalysisStruct(:).rms] ; 

p_K_i = kruskalwallis(K_i,flyLine,'off') ;
p_K_p = kruskalwallis(K_p,flyLine,'off') ;
p_gainRatio = kruskalwallis(gainRatio,flyLine,'off') ;
p_deltaT = kruskalwallis(deltaT,flyLine,'off') ;
p_rms = kruskalwallis(deltaT,flyLine,'off') ;

disp('p value for K_i:')
disp(p_K_i)
disp('p value for K_p:')
disp(p_K_p)
disp('p value for gain ratio:')
disp(p_gainRatio)
disp('p value for \Delta T:')
disp(p_deltaT)
disp('p value for rms:')
disp(p_rms)



