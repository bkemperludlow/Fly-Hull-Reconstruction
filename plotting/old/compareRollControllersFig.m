cd('G:\Janelia Flies\kir2.1 flies round 2\Analysis\Roll Controller Analysis')
load controllerAnalysisStruct

%find indices for relevant things to plot
controlRightInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == 2) ;
controlLeftInd = find([controllerAnalysisStruct(:).flyType] == 2 & ...
    [controllerAnalysisStruct(:).pertType] == -2) ;
experimentalRightInd = find([controllerAnalysisStruct(:).flyType] == 1 & ...
    [controllerAnalysisStruct(:).pertType] == 2) ;
experimentalLeftInd = find([controllerAnalysisStruct(:).flyType] == 1 & ...
    [controllerAnalysisStruct(:).pertType] == -2) ;

controlInd = find([controllerAnalysisStruct(:).flyType] == 2) ;
experimentalInd = find([controllerAnalysisStruct(:).flyType] == 1) ;

%make some plots:
%=========================================================================
%gain ratio
h_gain = figure('Position',[680   302   193  676]) ; 
set(gcf,'PaperPositionMode','auto')
%-------------------------------------------------------
subplot(3,1,1)
ylim1 = [-0.01+ min([controllerAnalysisStruct(:).K_i] ./ ...
    (1000*[controllerAnalysisStruct(:).K_p])), 0.01 + max([controllerAnalysisStruct(:).K_i] ./ ...
    (1000*[controllerAnalysisStruct(:).K_p]))] ; 
%ylim1 = [-0.12 0.2] ; 
hold on

plot(2*ones(size(controlRightInd)) + 0.1*(rand(size(controlRightInd)) - 0.5) + 0.1,...
    [controllerAnalysisStruct(controlRightInd).K_i] ./ ...
    (1000*[controllerAnalysisStruct(controlRightInd).K_p]), 'k^','MarkerFaceColor',[0 .7 0])
plot(ones(size(experimentalRightInd))+ 0.1*(rand(size(experimentalRightInd)) - 0.5) + 0.1, ...
    [controllerAnalysisStruct(experimentalRightInd).K_i] ./ ...
    (1000*[controllerAnalysisStruct(experimentalRightInd).K_p]), 'k^','MarkerFaceColor',[.7 0 0])

plot(2*ones(size(controlLeftInd)) + 0.1*(rand(size(controlLeftInd)) - 0.5) - 0.1, ...
    [controllerAnalysisStruct(controlLeftInd).K_i] ./ ...
    (1000*[controllerAnalysisStruct(controlLeftInd).K_p]), 'kv','MarkerFaceColor',[0 .7 0])
plot(ones(size(experimentalLeftInd)) + 0.1*(rand(size(experimentalLeftInd)) - 0.5) - 0.1, ...
    [controllerAnalysisStruct(experimentalLeftInd).K_i] ./ ...
    (1000*[controllerAnalysisStruct(experimentalLeftInd).K_p]), 'kv','MarkerFaceColor',[.7 0 0])

set(gca,'xlim',[.5 2.5]) 
set(gca,'ylim',ylim1) 
ylabel('K_i / K_p [1/ms]')

%-------------------------------------------------------------------------

subplot(3,1,2)
ylim2 = [-0.1+ min([controllerAnalysisStruct(:).K_i]), 0.1 + max([controllerAnalysisStruct(:).K_i])] ; 
hold on

plot(2*ones(size(controlRightInd)) + 0.1*(rand(size(controlRightInd)) - 0.5) + 0.1, ...
    [controllerAnalysisStruct(controlRightInd).K_i],...
    'k^','MarkerFaceColor',[0 .7 0])
plot(ones(size(experimentalRightInd)) + 0.1*(rand(size(experimentalRightInd)) - 0.5) + 0.1, ...
    [controllerAnalysisStruct(experimentalRightInd).K_i], ...
    'k^','MarkerFaceColor',[.7 0 0])

plot(2*ones(size(controlLeftInd)) + 0.1*(rand(size(controlLeftInd)) - 0.5) - 0.1,...
    [controllerAnalysisStruct(controlLeftInd).K_i],...
    'kv','MarkerFaceColor',[0 .7 0])
plot(ones(size(experimentalLeftInd)) + 0.1*(rand(size(experimentalLeftInd)) - 0.5) - 0.1, ...
    [controllerAnalysisStruct(experimentalLeftInd).K_i],...
    'kv','MarkerFaceColor',[.7 0 0])
set(gca,'xlim',[.5 2.5]) 
set(gca,'ylim',ylim2)
ylabel('K_i [deg]')

%-------------------------------------------------------------------------

subplot(3,1,3)
ylim3 = [-0.1+ min(1000*[controllerAnalysisStruct(:).K_p]), 0.1 + max(1000*[controllerAnalysisStruct(:).K_p])] ; 

hold on
plot(2*ones(size(controlRightInd)) + 0.1*(rand(size(controlRightInd)) - 0.5) + 0.1,...
    1000*[controllerAnalysisStruct(controlRightInd).K_p],...
    'k^','MarkerFaceColor',[0 .7 0])
plot(ones(size(experimentalRightInd)) + 0.1*(rand(size(experimentalRightInd)) - 0.5) + 0.1, ...
    1000*[controllerAnalysisStruct(experimentalRightInd).K_p], ...
    'k^','MarkerFaceColor',[.7 0 0])

plot(2*ones(size(controlLeftInd)) + 0.1*(rand(size(controlLeftInd)) - 0.5) - 0.1, ...
    1000*[controllerAnalysisStruct(controlLeftInd).K_p],...
    'kv','MarkerFaceColor',[0 .7 0])
plot(ones(size(experimentalLeftInd)) + 0.1*(rand(size(experimentalLeftInd)) - 0.5) - 0.1, ...
    1000*[controllerAnalysisStruct(experimentalLeftInd).K_p],...
    'kv','MarkerFaceColor',[.7 0 0])
set(gca,'xlim',[.5 2.5]) 
set(gca,'ylim',ylim3)
ylabel('K_p [ms]')

%=========================================================================
%scatter

h_scatter = figure ; 
set(gcf,'PaperPositionMode','auto')
hold on

plot([controllerAnalysisStruct(controlRightInd).P_norm], ...
    [controllerAnalysisStruct(controlRightInd).I_norm], 'k^', 'MarkerFaceColor', [0 .7 0])
plot([controllerAnalysisStruct(controlLeftInd).P_norm], ...
    [controllerAnalysisStruct(controlLeftInd).I_norm], 'kv', 'MarkerFaceColor', [0 .7 0])

plot([controllerAnalysisStruct(experimentalRightInd).P_norm], ...
    [controllerAnalysisStruct(experimentalRightInd).I_norm], 'k^', 'MarkerFaceColor', [.7 0 0])
plot([controllerAnalysisStruct(experimentalLeftInd).P_norm], ...
    [controllerAnalysisStruct(experimentalLeftInd).I_norm], 'kv', 'MarkerFaceColor', [.7 0 0])
ylabel('I norm')
xlabel('P norm')
axis tight

h_scatter2 = figure ; 
set(gcf,'PaperPositionMode','auto')
hold on

plot([controllerAnalysisStruct(controlRightInd).K_p], ...
    [controllerAnalysisStruct(controlRightInd).K_i], 'k^', 'MarkerFaceColor', [0 .7 0])
plot([controllerAnalysisStruct(controlLeftInd).K_p], ...
    [controllerAnalysisStruct(controlLeftInd).K_i], 'kv', 'MarkerFaceColor', [0 .7 0])

plot([controllerAnalysisStruct(experimentalRightInd).K_p], ...
    [controllerAnalysisStruct(experimentalRightInd).K_i], 'k^', 'MarkerFaceColor', [.7 0 0])
plot([controllerAnalysisStruct(experimentalLeftInd).K_p], ...
    [controllerAnalysisStruct(experimentalLeftInd).K_i], 'kv', 'MarkerFaceColor', [.7 0 0])
ylabel('K_i')
xlabel('K_p')
axis tight

%=========================================================================
% histogram
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

%{
cd('G:\b1 paper\michael meeting\')
print(h_gain,'controllerContributions_gain.svg','-dsvg')
print(h_bar,'controllerContributions_hist_gain.svg','-dsvg')
print(h_deltaT,'deltaT_hist.svg','-dsvg')


%}