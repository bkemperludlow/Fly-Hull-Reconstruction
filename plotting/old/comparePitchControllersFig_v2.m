cd('G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Controller Analysis')
load controllerAnalysisStruct

%find indices for relevant things to plot
MB204B_w1118_UpInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == 1 & ([controllerAnalysisStruct(:).ExprNum] == 2 | ...
    [controllerAnalysisStruct(:).ExprNum] == 7)) ;
MB204B_w1118_DownInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == -1 & ([controllerAnalysisStruct(:).ExprNum] == 2 | ...
    [controllerAnalysisStruct(:).ExprNum] == 7)) ;
MB204B_kir_UpInd = find([controllerAnalysisStruct(:).flyType] == 1 & ...
    [controllerAnalysisStruct(:).pertType] == 1 & ([controllerAnalysisStruct(:).ExprNum] == 1 | ...
    [controllerAnalysisStruct(:).ExprNum] == 8)) ;
MB204B_kir_DownInd = find([controllerAnalysisStruct(:).flyType] == 1 & ...
    [controllerAnalysisStruct(:).pertType] == -1 & ([controllerAnalysisStruct(:).ExprNum] == 1 | ...
    [controllerAnalysisStruct(:).ExprNum] == 8)) ;

EmptyVector_UpInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == 1 & ([controllerAnalysisStruct(:).ExprNum] == 6 | ...
    [controllerAnalysisStruct(:).ExprNum] == 14 | [controllerAnalysisStruct(:).ExprNum] == 15 )) ;
EmptyVector_DownInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == -1 & ([controllerAnalysisStruct(:).ExprNum] == 6 | ...
    [controllerAnalysisStruct(:).ExprNum] == 14 | [controllerAnalysisStruct(:).ExprNum] == 15 )) ;

w1118_UpInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == 1 &  [controllerAnalysisStruct(:).ExprNum] == 5) ;
w1118_DownInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == -1 &  [controllerAnalysisStruct(:).ExprNum] == 5) ;

controlInd = find([controllerAnalysisStruct(:).flyType] == 2) ;
experimentalInd = find([controllerAnalysisStruct(:).flyType] == 1) ;

MB204B_w1118_Color = [0 .7 0 ] ; 
MB204B_kir_Color = [0.7 0 0 ] ; 
EmptyVector_Color = [238,232,170]/255 ; 
w1118_Color = [173,216,230]/255 ; 

IndCell = cell(4,2) ; 
IndCell{1,1} = MB204B_kir_UpInd ; IndCell{1,2} = MB204B_kir_DownInd ;
IndCell{2,1} = MB204B_w1118_UpInd ; IndCell{2,2} = MB204B_w1118_DownInd ;
IndCell{3,1} = EmptyVector_UpInd ; IndCell{3,2} = EmptyVector_DownInd ;
IndCell{4,1} = w1118_UpInd ; IndCell{4,2} = w1118_DownInd ;

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

%pitch up controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,1} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) + 0.1,...
        [controllerAnalysisStruct(indTemp).K_i] ./ ...
        (1000*[controllerAnalysisStruct(indTemp).K_p]), 'k^','MarkerFaceColor',colorTemp)
end

%pitch down controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,2} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) - 0.1,...
        [controllerAnalysisStruct(indTemp).K_i] ./ ...
        (1000*[controllerAnalysisStruct(indTemp).K_p]), 'kv','MarkerFaceColor',colorTemp)
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

%pitch up controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,1} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) + 0.1,...
        [controllerAnalysisStruct(indTemp).K_i], 'k^','MarkerFaceColor',colorTemp)
end

%pitch down controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,2} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) - 0.1,...
        [controllerAnalysisStruct(indTemp).K_i], 'kv','MarkerFaceColor',colorTemp)
end
set(gca,'xlim',xlim) 
set(gca,'ylim',ylim2)
ylabel('K_i [deg]')

%-------------------------------------------------------------------------

subplot(2,2,3)
%ylim3 = [-0.1+ min(1000*[controllerAnalysisStruct(:).K_p]), 0.1 + max(1000*[controllerAnalysisStruct(:).K_p])] ; 
ylim3 = [0, 0.1 + max(1000*[controllerAnalysisStruct(:).K_p])] ; 

hold on
%pitch up controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,1} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) + 0.1,...
        1000*[controllerAnalysisStruct(indTemp).K_p], 'k^','MarkerFaceColor',colorTemp)
end

%pitch down controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,2} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) - 0.1,...
        1000*[controllerAnalysisStruct(indTemp).K_p], 'kv','MarkerFaceColor',colorTemp)
end
set(gca,'xlim',xlim) 
set(gca,'ylim',ylim3)
ylabel('K_p [ms]')

%--------------------------------------------------------------------------------------------

subplot(2,2,4)
%ylim3 = [-0.1+ min(1000*[controllerAnalysisStruct(:).K_p]), 0.1 + max(1000*[controllerAnalysisStruct(:).K_p])] ; 
ylim4 = [-0.1 + min(1000*[controllerAnalysisStruct(:).deltaT]), 0.1 + max(1000*[controllerAnalysisStruct(:).deltaT])] ; 

hold on
%pitch up controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,1} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) + 0.1,...
        1000*[controllerAnalysisStruct(indTemp).deltaT], 'k^','MarkerFaceColor',colorTemp)
end

%pitch down controllers
for i = 1:size(IndCell,1)
    indTemp = IndCell{i,2} ; 
    colorTemp = ColorMat(i, :) ; 
    plot(i*ones(size(indTemp)) + 0.1*(rand(size(indTemp)) - 0.5) - 0.1,...
        1000*[controllerAnalysisStruct(indTemp).deltaT], 'kv','MarkerFaceColor',colorTemp)
end
set(gca,'xlim',xlim) 
set(gca,'ylim',ylim4)
ylabel('\Delta T [ms]')

%=========================================================================
%scatter
%{
h_scatter = figure ; 
set(gcf,'PaperPositionMode','auto')
hold on

plot([controllerAnalysisStruct(MB204B_w1118_UpInd).P_norm], ...
    [controllerAnalysisStruct(MB204B_w1118_UpInd).I_norm], 'k^', 'MarkerFaceColor', [0 .7 0])
plot([controllerAnalysisStruct(MB204B_w1118_DownInd).P_norm], ...
    [controllerAnalysisStruct(MB204B_w1118_DownInd).I_norm], 'kv', 'MarkerFaceColor', [0 .7 0])

plot([controllerAnalysisStruct(MB204B_kir_UpInd).P_norm], ...
    [controllerAnalysisStruct(MB204B_kir_UpInd).I_norm], 'k^', 'MarkerFaceColor', [.7 0 0])
plot([controllerAnalysisStruct(MB204B_kir_DownInd).P_norm], ...
    [controllerAnalysisStruct(MB204B_kir_DownInd).I_norm], 'kv', 'MarkerFaceColor', [.7 0 0])
ylabel('I norm')
xlabel('P norm')
axis tight

h_scatter2 = figure ; 
set(gcf,'PaperPositionMode','auto')
hold on

plot([controllerAnalysisStruct(MB204B_w1118_UpInd).K_p], ...
    [controllerAnalysisStruct(MB204B_w1118_UpInd).K_i], 'k^', 'MarkerFaceColor', [0 .7 0])
plot([controllerAnalysisStruct(MB204B_w1118_DownInd).K_p], ...
    [controllerAnalysisStruct(MB204B_w1118_DownInd).K_i], 'kv', 'MarkerFaceColor', [0 .7 0])

plot([controllerAnalysisStruct(MB204B_kir_UpInd).K_p], ...
    [controllerAnalysisStruct(MB204B_kir_UpInd).K_i], 'k^', 'MarkerFaceColor', [.7 0 0])
plot([controllerAnalysisStruct(MB204B_kir_DownInd).K_p], ...
    [controllerAnalysisStruct(MB204B_kir_DownInd).K_i], 'kv', 'MarkerFaceColor', [.7 0 0])
ylabel('K_i')
xlabel('K_p')
axis tight
%}
%=========================================================================
% histogram
%{
h_bar = figure('Position',[680   302   193*.6  676]) ; 
set(gcf,'PaperPositionMode','auto')
nbins = 3 ;
%----------------------------------------------------------------------
[f1_control, x1_control] = hist([controllerAnalysisStruct(controlInd).K_i] ./ ...
    (1000*[controllerAnalysisStruct(controlInd).K_p]), linspace(ylim1(1),ylim1(2),nbins)) ; 
[f1_experimental, x1_experimental] = hist([controllerAnalysisStruct(experimentalInd).K_i] ./ ...
    (1000*[controllerAnalysisStruct(experimentalInd).K_p]), linspace(ylim1(1),ylim1(2),nbins)) ; 

subplot(3,2,1)
hold on
barh( x1_control , f1_control, 'FaceColor',[0 .7 0 ]) 
set(gca,'ylim',ylim1)
subplot(3,2,2)
hold on
barh( x1_experimental , f1_experimental, 'FaceColor',[.7 0 0 ]) 
set(gca,'ylim',ylim1)
%----------------------------------------------------------------------
[f2_control, x2_control] = hist([controllerAnalysisStruct(controlInd).K_i],...
    linspace(ylim2(1),ylim2(2),nbins)) ; 
[f2_experimental, x2_experimental] = hist([controllerAnalysisStruct(experimentalInd).K_i],...
    linspace(ylim2(1),ylim2(2),nbins)) ; 

subplot(3,2,3)
hold on
barh( x2_control , f2_control, 'FaceColor',[0 .7 0 ]) 
set(gca,'ylim',ylim2)
subplot(3,2,4)
hold on
barh( x2_experimental , f2_experimental, 'FaceColor',[.7 0 0 ]) 
set(gca,'ylim',ylim2)
%----------------------------------------------------------------------
[f3_control, x3_control] = hist(1000*[controllerAnalysisStruct(controlInd).K_p],...
    linspace(ylim3(1),ylim3(2),nbins)) ; 
[f3_experimental, x3_experimental] = hist(1000*[controllerAnalysisStruct(experimentalInd).K_p],...
    linspace(ylim3(1),ylim3(2),nbins)) ; 

subplot(3,2,5)
hold on
barh( x3_control , f3_control, 'FaceColor',[0 .7 0 ]) 
set(gca,'ylim',ylim3)
subplot(3,2,6)
hold on
barh( x3_experimental , f3_experimental, 'FaceColor',[.7 0 0 ]) 
set(gca,'ylim',ylim3)
%}

%=========================================================================
%delta t distrtibution
%{
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