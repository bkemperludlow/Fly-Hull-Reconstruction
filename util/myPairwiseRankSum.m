cd('G:\Janelia Flies\kir2.1 flies round 2\Analysis\Roll Controller Analysis')
load controllerAnalysisStruct

controllerVariable = 'deltaT' ; %'gainRatio', 'K_i', 'K_p', or 'deltaT'

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

labelList = {'MB204B > kir2.1', 'MB204B > w1118', 'EmptyVector > kir2.1', 'w1118 > kir2.1'} ;

IndCell = cell(4,1) ; 
IndCell{1,1} = MB204B_kir_Ind ; 
IndCell{2,1} = MB204B_w1118_Ind ; 
IndCell{3,1} = EmptyVector_Ind ; 
IndCell{4,1} = w1118_Ind ; 


K_i = [controllerAnalysisStruct(:).K_i] ;
K_p = [controllerAnalysisStruct(:).K_p] ;
gainRatio = K_i ./ K_p ;
deltaT = [controllerAnalysisStruct(:).deltaT] ;

switch controllerVariable
    case 'gainRatio'
        compareVar = gainRatio ; 
    case 'K_i' 
        compareVar = K_i ; 
    case 'K_p' 
        compareVar = K_p ; 
    case 'deltaT'
        compareVar = deltaT ; 
    otherwise
        keyboard ; 
end

pValueMat = nan(4) ; 
for i = 1:4
    for j = 1:4
        if i ==j 
            continue
        end
        
        temp1 = compareVar(IndCell{i}) ; 
        temp2 = compareVar(IndCell{j}) ; 
        pValueTemp =  ranksum(temp1, temp2) ; 
        pValueMat(i,j) = pValueTemp ; 
        
    end
end

cmax = max(abs(pValueMat(:) - 0.05)) ; 
figure ; imagesc(pValueMat - 0.05) 
colormap(redblue)
caxis([-cmax, cmax])
colorbar
set(gca,'XTick',1:4,'XTickLabel', labelList)
set(gca,'YTick',1:4,'YTickLabel', labelList)
set(gca,'XTickLabelRotation',45)
        
