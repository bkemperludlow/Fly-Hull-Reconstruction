cd('G:\Janelia Flies\kir2.1 flies round 2\Analysis\Roll Controller Analysis')
load accelerationAndPhiAmpStruct 

%% find indices for relevant things to plot
MB204B_w1118_RightInd = find([accelerationAndPhiAmpStruct(:).flyType] == 2 & ...
    [accelerationAndPhiAmpStruct(:).pertType] == 2 & ([accelerationAndPhiAmpStruct(:).ExprNum] == 2 | ...
    [accelerationAndPhiAmpStruct(:).ExprNum] == 7)) ;
MB204B_w1118_LeftInd = find([accelerationAndPhiAmpStruct(:).flyType] == 2 & ...
    [accelerationAndPhiAmpStruct(:).pertType] == -2 & ([accelerationAndPhiAmpStruct(:).ExprNum] == 2 | ...
    [accelerationAndPhiAmpStruct(:).ExprNum] == 7)) ;
MB204B_kir_RightInd = find([accelerationAndPhiAmpStruct(:).flyType] == 1 & ...
    [accelerationAndPhiAmpStruct(:).pertType] == 2 & ([accelerationAndPhiAmpStruct(:).ExprNum] == 1 | ...
    [accelerationAndPhiAmpStruct(:).ExprNum] == 8)) ;
MB204B_kir_LeftInd = find([accelerationAndPhiAmpStruct(:).flyType] == 1 & ...
    [accelerationAndPhiAmpStruct(:).pertType] == -2 & ([accelerationAndPhiAmpStruct(:).ExprNum] == 1 | ...
    [accelerationAndPhiAmpStruct(:).ExprNum] == 8)) ;

EmptyVector_RightInd = find([accelerationAndPhiAmpStruct(:).flyType] == 2 & ...
    [accelerationAndPhiAmpStruct(:).pertType] == 2 & ([accelerationAndPhiAmpStruct(:).ExprNum] == 6 | ...
    [accelerationAndPhiAmpStruct(:).ExprNum] == 14 | [accelerationAndPhiAmpStruct(:).ExprNum] == 15 )) ;
EmptyVector_LeftInd = find([accelerationAndPhiAmpStruct(:).flyType] == 2 & ...
    [accelerationAndPhiAmpStruct(:).pertType] == -2 & ([accelerationAndPhiAmpStruct(:).ExprNum] == 6 | ...
    [accelerationAndPhiAmpStruct(:).ExprNum] == 14 | [accelerationAndPhiAmpStruct(:).ExprNum] == 15 )) ;

w1118_RightInd = find([accelerationAndPhiAmpStruct(:).flyType] == 2 & ...
    [accelerationAndPhiAmpStruct(:).pertType] == 2 &  [accelerationAndPhiAmpStruct(:).ExprNum] == 5) ;
w1118_LeftInd = find([accelerationAndPhiAmpStruct(:).flyType] == 2 & ...
    [accelerationAndPhiAmpStruct(:).pertType] == -2 &  [accelerationAndPhiAmpStruct(:).ExprNum] == 5) ;

controlInd = find([accelerationAndPhiAmpStruct(:).flyType] == 2) ;
experimentalInd = find([accelerationAndPhiAmpStruct(:).flyType] == 1) ;

MB204B_w1118_Color = [0 .7 0 ] ; 
MB204B_kir_Color = [0.7 0 0 ] ; 
EmptyVector_Color = [238,232,170]/255 ; 
w1118_Color = [173,216,230]/255 ; 

IndCell = cell(4,2) ; 
IndCell{1,1} = MB204B_kir_RightInd ; IndCell{1,2} = MB204B_kir_LeftInd ;
IndCell{2,1} = MB204B_w1118_RightInd ; IndCell{2,2} = MB204B_w1118_LeftInd ;
IndCell{3,1} = EmptyVector_RightInd ; IndCell{3,2} = EmptyVector_LeftInd ;
IndCell{4,1} = w1118_RightInd ; IndCell{4,2} = w1118_LeftInd ;

ColorMat = [MB204B_kir_Color ; MB204B_w1118_Color ; EmptyVector_Color ; w1118_Color] ;  

labelList = {'MB204B > kir2.1', 'MB204B > w1118', 'EmptyVector > kir2.1', 'w1118 > kir2.1'} ;


badInd = find(([accelerationAndPhiAmpStruct(:).ExprNum] == 1 & ...
    [accelerationAndPhiAmpStruct(:).MovNum] == 65) | ([accelerationAndPhiAmpStruct(:).ExprNum] == 5 & ...
    [accelerationAndPhiAmpStruct(:).MovNum] == 1)) ; 

maxRollAccel = [accelerationAndPhiAmpStruct(:).maxRollAccel] ; 
maxPhiAmpDiff = [accelerationAndPhiAmpStruct(:).maxPhiAmpDiff] ; 

maxRollAccel(badInd) = nan ;
deltaPhiFront(badInd) = nan ;
deltaPhiBack(badInd) = nan ;

%% make figure
h_main = figure ; 
hold on

%roll right controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,1} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(maxPhiAmpDiff(indTemp), maxRollAccel(indTemp), 'k>','MarkerFaceColor',colorTemp)
end

%roll left controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,2} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(maxPhiAmpDiff(indTemp), maxRollAccel(indTemp), 'k<','MarkerFaceColor',colorTemp)
end


xlabel('\Delta\Phi_{diff} [deg]')
ylabel('Roll Accel. [deg/s^2]')
axis tight
nanInd = isnan(maxRollAccel) |  isnan(maxPhiAmpDiff) ;
R = corrcoef(maxRollAccel(~nanInd),...
    maxPhiAmpDiff(~nanInd)') ; 
R_sq_front = R(1,2)^2 ;
title(['R^2 = ' num2str(R_sq_front)]) 