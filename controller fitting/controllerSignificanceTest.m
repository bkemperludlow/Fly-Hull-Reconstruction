function [significanceStruct, pValueLabels] = controllerSignificanceTest()
% want to have p-values to report for the differences between b1-silenced
% and control group flies

cd('G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Controller Analysis')
load controllerAnalysisStruct

significanceStruct = struct('ComparisonType',[], 'pValues_tt', [], 'pValues_rs',[],...
    'pValues_ks',[], 'N1', [], 'N2', []) ;
pValueLabels = {'Norm Ratio','I Norm', 'P Norm', 'Gain Ratio', 'K_i', 'K_p', 'deltaT'} ;

cc = 1 ;

controlUpInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pitchType] == 1) ;
controlDownInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pitchType] == -1) ;
experimentalUpInd = find([controllerAnalysisStruct(:).flyType] == 1 & ...
    [controllerAnalysisStruct(:).pitchType] == 1) ;
experimentalDownInd = find([controllerAnalysisStruct(:).flyType] == 1 & ...
    [controllerAnalysisStruct(:).pitchType] == -1) ;

controlInd = find([controllerAnalysisStruct(:).flyType] == 2) ;
experimentalInd = find([controllerAnalysisStruct(:).flyType] == 1) ;


%% all control vs all experimental 1

Ind1 = controlInd ; 
Ind2 = experimentalInd ; 

p_all_rs = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'rs') ;
p_all_tt = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'tt') ;
p_all_ks = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'ks') ;

significanceStruct(cc).ComparisonType = 'All control vs all experimental' ;
significanceStruct(cc).pValues_tt = p_all_tt ;
significanceStruct(cc).pValues_rs = p_all_rs ;
significanceStruct(cc).pValues_ks = p_all_ks ;
significanceStruct(cc).N1 = length(Ind1) ;
significanceStruct(cc).N2 = length(Ind2) ;
cc = cc + 1 ; 

%% control pitch up vs down 2

Ind1 = controlUpInd ; 
Ind2 = controlDownInd ; 

p_all_rs = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'rs') ;
p_all_tt = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'tt') ;
p_all_ks = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'ks') ;

significanceStruct(cc).ComparisonType = 'Control pitch up vs control pitch down' ;
significanceStruct(cc).pValues_tt = p_all_tt ;
significanceStruct(cc).pValues_rs = p_all_rs ;
significanceStruct(cc).pValues_ks = p_all_ks ;
significanceStruct(cc).N1 = length(Ind1) ;
significanceStruct(cc).N2 = length(Ind2) ;
cc = cc + 1 ; 

%% experimental pitch up vs down 3
Ind1 = experimentalUpInd ; 
Ind2 = experimentalDownInd ; 

p_all_rs = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'rs') ;
p_all_tt = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'tt') ;
p_all_ks = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'ks') ;

significanceStruct(cc).ComparisonType = 'Experimental pitch up vs experimental pitch down' ;
significanceStruct(cc).pValues_tt = p_all_tt ;
significanceStruct(cc).pValues_rs = p_all_rs ;
significanceStruct(cc).pValues_ks = p_all_ks ;
significanceStruct(cc).N1 = length(Ind1) ;
significanceStruct(cc).N2 = length(Ind2) ;
cc = cc + 1 ; 

%% control pitch down vs experimental pitch down 4
Ind1 = controlDownInd ; 
Ind2 = experimentalDownInd ; 

p_all_rs = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'rs') ;
p_all_tt = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'tt') ;
p_all_ks = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'ks') ;

significanceStruct(cc).ComparisonType = 'Control pitch down vs experimental pitch down' ;
significanceStruct(cc).pValues_tt = p_all_tt ;
significanceStruct(cc).pValues_rs = p_all_rs ;
significanceStruct(cc).pValues_ks = p_all_ks ;
significanceStruct(cc).N1 = length(Ind1) ;
significanceStruct(cc).N2 = length(Ind2) ;
cc = cc + 1 ; 

%% control pitch up vs experimental pitch up 5
Ind1 = controlUpInd ; 
Ind2 = experimentalUpInd ; 

p_all_rs = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'rs') ;
p_all_tt = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'tt') ;
p_all_ks = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'ks') ;

significanceStruct(cc).ComparisonType = 'Control pitch up vs experimental pitch up' ;
significanceStruct(cc).pValues_tt = p_all_tt ;
significanceStruct(cc).pValues_rs = p_all_rs ;
significanceStruct(cc).pValues_ks = p_all_ks ;
significanceStruct(cc).N1 = length(Ind1) ;
significanceStruct(cc).N2 = length(Ind2) ;
cc = cc + 1 ; 

%{
%% MB258C (both up and down) 6
Ind1 = find([controllerAnalysisStruct(:).ExprNum] == 7 | [controllerAnalysisStruct(:).ExprNum] == 12 ) ;
Ind2 = find([controllerAnalysisStruct(:).ExprNum] == 6 | [controllerAnalysisStruct(:).ExprNum] == 13 ) ;

p_all_rs = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'rs') ;
p_all_tt = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'tt') ;
p_all_ks = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'ks') ;

significanceStruct(cc).ComparisonType = 'Control vs experimental for MB258C (both up and down)' ;
significanceStruct(cc).pValues_tt = p_all_tt ;
significanceStruct(cc).pValues_rs = p_all_rs ;
significanceStruct(cc).pValues_ks = p_all_ks ;
significanceStruct(cc).N1 = length(Ind1) ;
significanceStruct(cc).N2 = length(Ind2) ;
cc = cc + 1 ; 

%% MB204B (both up and down) 7
Ind1 = find([controllerAnalysisStruct(:).ExprNum] == 4 ) ;
Ind2 = find([controllerAnalysisStruct(:).ExprNum] == 3 | [controllerAnalysisStruct(:).ExprNum] == 8 ) ;

p_all_rs = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'rs') ;
p_all_tt = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'tt') ;
p_all_ks = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'ks') ;

significanceStruct(cc).ComparisonType = 'Control vs experimental for MB204B (both up and down)'  ;
significanceStruct(cc).pValues_tt = p_all_tt ;
significanceStruct(cc).pValues_rs = p_all_rs ;
significanceStruct(cc).pValues_ks = p_all_ks ;
significanceStruct(cc).N1 = length(Ind1) ;
significanceStruct(cc).N2 = length(Ind2) ;
cc = cc + 1 ; 

%% 9281(2) (both up and down) 8
Ind1 = find([controllerAnalysisStruct(:).ExprNum] == 11 ) ;
Ind2 = find([controllerAnalysisStruct(:).ExprNum] == 9 | [controllerAnalysisStruct(:).ExprNum] == 10 ) ;

p_all_rs = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'rs') ;
p_all_tt = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'tt') ;
p_all_ks = panelSignificanceTest(controllerAnalysisStruct, Ind1, Ind2, 'ks') ;

significanceStruct(cc).ComparisonType = 'Control vs experimental for 9281(2) (both up and down)'  ;
significanceStruct(cc).pValues_tt = p_all_tt ;
significanceStruct(cc).pValues_rs = p_all_rs ;
significanceStruct(cc).pValues_ks = p_all_ks ;
significanceStruct(cc).N1 = length(Ind1) ;
significanceStruct(cc).N2 = length(Ind2) ;
cc = cc + 1 ; 
%}

end

%% do the panel of tests given a set of indices
function pValues = panelSignificanceTest(controllerAnalysisStruct, ind1, ind2, testType)

NormRatio1 = [controllerAnalysisStruct(ind1).I_norm] ./ ...
    [controllerAnalysisStruct(ind1).P_norm] ; 
NormRatio2 = [controllerAnalysisStruct(ind2).I_norm] ./ ...
    [controllerAnalysisStruct(ind2).P_norm] ; 

INorm1 = [controllerAnalysisStruct(ind1).I_norm]  ; 
INorm2 = [controllerAnalysisStruct(ind2).I_norm]  ;

PNorm1 = [controllerAnalysisStruct(ind1).P_norm]  ; 
PNorm2 = [controllerAnalysisStruct(ind2).P_norm]  ;

GainRatio1 = [controllerAnalysisStruct(ind1).K_i] ./ ...
    [controllerAnalysisStruct(ind1).K_p]  ; 
GainRatio2 = [controllerAnalysisStruct(ind2).K_i]  ./ ...
    [controllerAnalysisStruct(ind2).K_p]  ;

K_i1 = [controllerAnalysisStruct(ind1).K_i]  ; 
K_i2 = [controllerAnalysisStruct(ind2).K_i]  ; 

K_p1 = [controllerAnalysisStruct(ind1).K_p]  ; 
K_p2 = [controllerAnalysisStruct(ind2).K_p]  ; 

deltaT1 = ([controllerAnalysisStruct(ind1).deltaT]./ ...
    [controllerAnalysisStruct(ind1).medianWingBeat]) ;
deltaT2 = ([controllerAnalysisStruct(ind2).deltaT]./ ...
    [controllerAnalysisStruct(ind2).medianWingBeat]) ;


switch testType
    %rank sum test
    case 'rs'
        [p_NormRatio, ~] = ranksum(NormRatio1, NormRatio2) ;
        [p_INorm, ~] = ranksum(INorm1, INorm2) ;
        [p_PNorm, ~] = ranksum(PNorm1, PNorm2) ;
        [p_GainRatio, ~] = ranksum(GainRatio1, GainRatio2) ;
        [p_K_i, ~] = ranksum(K_i1, K_i2) ;
        [p_K_p, ~] = ranksum(K_p1, K_p2) ;
        [p_deltaT, ~] = ranksum(deltaT1, deltaT2) ;
    % two sample t test
    case 'tt'
        [~, p_NormRatio] = ttest2(NormRatio1, NormRatio2) ;
        [~, p_INorm] = ttest2(INorm1, INorm2) ;
        [~, p_PNorm] = ttest2(PNorm1, PNorm2) ;
        [~, p_GainRatio] = ttest2(GainRatio1, GainRatio2) ;
        [~, p_K_i] = ttest2(K_i1, K_i2) ;
        [~, p_K_p] = ttest2(K_p1, K_p2) ;
        [~, p_deltaT] = ttest2(deltaT1, deltaT2) ;
    % kolmogorov-smirnov test
    case 'ks'
        [~, p_NormRatio] = kstest2(NormRatio1, NormRatio2) ;
        [~, p_INorm] = kstest2(INorm1, INorm2) ;
        [~, p_PNorm] = kstest2(PNorm1, PNorm2) ;
        [~, p_GainRatio] = kstest2(GainRatio1, GainRatio2) ;
        [~, p_K_i] = kstest2(K_i1, K_i2) ;
        [~, p_K_p] = kstest2(K_p1, K_p2) ;
        [~, p_deltaT] = kstest2(deltaT1, deltaT2) ;
    otherwise
        p_NormRatio = nan ;
        p_INorm = nan ;
        p_PNorm = nan ;
        p_GainRatio = nan ;
        p_K_i = nan ;
        p_K_p = nan ;
        p_deltaT = nan ;
        keyboard ;
        
end

pValues = [p_NormRatio, p_INorm, p_PNorm, p_GainRatio, p_K_i, p_K_p, p_deltaT] ; 
end

