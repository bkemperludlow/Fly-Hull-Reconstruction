cd('G:\Janelia Flies\kir2.1 flies round 2\Analysis\Roll Controller Analysis')
load controllerAnalysisStruct

%find indices for relevant things to plot

MB204B_w1118_RightInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == 2 & ([controllerAnalysisStruct(:).ExprNum] == 2 | ...
    [controllerAnalysisStruct(:).ExprNum] == 7)) ;
MB204B_w1118_LeftInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == -2 & ([controllerAnalysisStruct(:).ExprNum] == 2 | ...
    [controllerAnalysisStruct(:).ExprNum] == 7)) ;
MB204B_kir_RightInd = find([controllerAnalysisStruct(:).flyType] == 1 & ...
    [controllerAnalysisStruct(:).pertType] == 2 & ([controllerAnalysisStruct(:).ExprNum] == 1 | ...
    [controllerAnalysisStruct(:).ExprNum] == 8)) ;
MB204B_kir_LeftInd = find([controllerAnalysisStruct(:).flyType] == 1 & ...
    [controllerAnalysisStruct(:).pertType] == -2 & ([controllerAnalysisStruct(:).ExprNum] == 1 | ...
    [controllerAnalysisStruct(:).ExprNum] == 8)) ;

EmptyVector_RightInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == 2 & ([controllerAnalysisStruct(:).ExprNum] == 6 | ...
    [controllerAnalysisStruct(:).ExprNum] == 14 | [controllerAnalysisStruct(:).ExprNum] == 15 )) ;
EmptyVector_LeftInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == -2 & ([controllerAnalysisStruct(:).ExprNum] == 6 | ...
    [controllerAnalysisStruct(:).ExprNum] == 14 | [controllerAnalysisStruct(:).ExprNum] == 15 )) ;

w1118_RightInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == 2 &  [controllerAnalysisStruct(:).ExprNum] == 5) ;
w1118_LeftInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == -2 &  [controllerAnalysisStruct(:).ExprNum] == 5) ;

controlInd = find([controllerAnalysisStruct(:).flyType] == 2) ;
experimentalInd = find([controllerAnalysisStruct(:).flyType] == 1) ;

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

%make some plots:
%=========================================================================
%gain ratio
h_gain = figure('Position',[680   302   2*193  (2/3)*676]) ; 
set(gcf,'PaperPositionMode','auto')
%-------------------------------------------------------
subplot(2,2,1)
ylim1 = [-0.01+ min([controllerAnalysisStruct(:).K_i] ./ ...
    (1000*[controllerAnalysisStruct(:).K_p])), 0.01 + max([controllerAnalysisStruct(:).K_i] ./ ...
    (1000*[controllerAnalysisStruct(:).K_p]))] ; 
xlim = [0.5 4.5] ; 
hold on

plot(linspace(xlim(1),xlim(2),50), zeros(1,50), '--', 'Color', .7*[1 1 1],'LineWidth',1.5)

%roll right controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,1} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) + 0.1,...
        [controllerAnalysisStruct(indTemp).K_i] ./ ...
        (1000*[controllerAnalysisStruct(indTemp).K_p]), 'k>','MarkerFaceColor',colorTemp)
end

%roll left controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,2} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) - 0.1,...
        [controllerAnalysisStruct(indTemp).K_i] ./ ...
        (1000*[controllerAnalysisStruct(indTemp).K_p]), 'k<','MarkerFaceColor',colorTemp)
end

set(gca,'xlim',xlim)
set(gca,'ylim',ylim1) 
ylabel('K_i / K_p [1/ms]')
%set(gca,'XTickLabel',labelList)
%set(gca,'XTickLabelRotation',45)
%-------------------------------------------------------------------------

subplot(2,2,2)
ylim2 = [-0.1+ min([controllerAnalysisStruct(:).K_i]), 0.1 + max([controllerAnalysisStruct(:).K_i])] ; 
hold on

plot(linspace(xlim(1),xlim(2),50), zeros(1,50), '--', 'Color', .7*[1 1 1],'LineWidth',1.5)

%roll right controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,1} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) + 0.1,...
        [controllerAnalysisStruct(indTemp).K_i], 'k>','MarkerFaceColor',colorTemp)
end

%roll left controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,2} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) - 0.1,...
        [controllerAnalysisStruct(indTemp).K_i], 'k<','MarkerFaceColor',colorTemp)
end
set(gca,'xlim',xlim) 
set(gca,'ylim',ylim2)
ylabel('K_i [deg]')

%-------------------------------------------------------------------------

subplot(2,2,3)
%ylim3 = [-0.1+ min(1000*[controllerAnalysisStruct(:).K_p]), 0.1 + max(1000*[controllerAnalysisStruct(:).K_p])] ; 
ylim3 = [0, 0.1 + max(1000*[controllerAnalysisStruct(:).K_p])] ; 

hold on
%roll right controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,1} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) + 0.1,...
        1000*[controllerAnalysisStruct(indTemp).K_p], 'k>','MarkerFaceColor',colorTemp)
end

%roll left controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,2} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) - 0.1,...
        1000*[controllerAnalysisStruct(indTemp).K_p], 'k<','MarkerFaceColor',colorTemp)
end
set(gca,'xlim',xlim) 
set(gca,'ylim',ylim3)
ylabel('K_p [ms]')

%--------------------------------------------------------------------------------------------

subplot(2,2,4)
%ylim3 = [-0.1+ min(1000*[controllerAnalysisStruct(:).K_p]), 0.1 + max(1000*[controllerAnalysisStruct(:).K_p])] ; 
ylim4 = [-0.1 + min(1000*[controllerAnalysisStruct(:).deltaT]), 0.1 + max(1000*[controllerAnalysisStruct(:).deltaT])] ; 

hold on
%roll right controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,1} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) + 0.1,...
        1000*[controllerAnalysisStruct(indTemp).deltaT], 'k>','MarkerFaceColor',colorTemp)
end

%roll left controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,2} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) - 0.1,...
        1000*[controllerAnalysisStruct(indTemp).deltaT], 'k<','MarkerFaceColor',colorTemp)
end
set(gca,'xlim',xlim) 
set(gca,'ylim',ylim4)
ylabel('\Delta T [ms]')

%{
%=========================================================================
%delta t distrtibution
nbins2 = 4 ;

deltaT_control = ([controllerAnalysisStruct(controlInd).deltaT]./ [controllerAnalysisStruct(controlInd).medianWingBeat]);
deltaT_experimental = ([controllerAnalysisStruct(experimentalInd).deltaT]./ [controllerAnalysisStruct(experimentalInd).medianWingBeat]);

[f, x] = hist([deltaT_control , deltaT_experimental], nbins2) ;
[f_ctrl, x_ctrl] = hist(deltaT_control,x) ;
[f_exp, x_exp] = hist(deltaT_experimental,x) ;

h_deltaT = figure ; 
set(gcf,'PaperPositionMode','auto')
subplot(2,1,1) 
bar(x_ctrl, f_ctrl, 'FaceColor',[0 .7 0]) 
ylabel('Counts')
subplot(2,1,2) 
bar(x_exp, f_exp, 'FaceColor',[.7 0 0]) 
xlabel('\Delta T [wingbeats]')
ylabel('Counts')
%}
%{
cd('G:\b1 paper\michael meeting\')
print(h_gain,'controllerContributions_gain.svg','-dsvg')
print(h_bar,'controllerContributions_hist_gain.svg','-dsvg')
print(h_deltaT,'deltaT_hist.svg','-dsvg')


%}