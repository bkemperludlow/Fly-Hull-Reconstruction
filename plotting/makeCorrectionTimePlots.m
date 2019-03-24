function [hresponse, hcorrection] = makeCorrectionTimePlots()

cd('F:\luca\Analysis\metadata\Correction Times')
correctionTimesUpData = importdata('CorrectionTimesUp11062014.mat') ;
correctionTimesDownData = importdata('CorrectionTimesDown11062014.mat') ;

t_l_up = 1000*correctionTimesUpData.T_latency ; %ms
t_l_down = 1000*correctionTimesDownData.T_latency ;

%{
arraySize = max([length(t_l_up),length(t_l_down)]) ; 
latencyWingbeatsUp = NaN(arraySize,1) ;
latencyWingbeatsDown = NaN(arraySize,1) ;

for i = 1:length(t_l_up)
    latencyWingbeatsUp(i) = (t_l_up(i) - mod(t_l_up(i),5))/5 ; 
end

for j = 1:length(t_l_down)
    latencyWingbeatsDown(j) = (t_l_down(j) - mod(t_l_down(j),5))/5 ; 
end
%}
%xbins = 0:3 ;
hresponse = figure('paperPositionMode','auto') ;
    set(gca,'fontsize',14) ;
    xlim = [4 13] ;
    ylim = [0 7] ;
    hold on
    
    avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
    tsfvec = [4.3 7.7 7.7 4.3 4.3] ; %mean + std for fit deltaT (6 +/- 1.7)
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 179 71]/255,'facealpha',1,'edgecolor','none') ;
    
    [f, x ] = hist([t_l_up;t_l_down;6.125;10.75],4) ; % movies ??? and 14.1 added
    %legend({'Pitch Up','Pitch Down'}, 'location','northwest')
    hbar = bar(x,f, 'edgecolor', 'k','facecolor',[0 .5 0] ) ;
    %{
    h = findobj(gca,'Type','patch');
    set(h(1),'FaceColor',[0 0 .8]);
    set(h(1),'EdgeColor','w');
    set(h(2),'FaceColor',[.8 0 0]);
    set(h(2),'EdgeColor','w');
    
    set(gca, 'Xtick', 0:3) ;
    set(gca, 'xlim',[0 3]) ; 
    %}
    set(gca,'xlim',xlim)
    xlabel('Correction Latency Time [ms]','fontsize',14)
    ylabel('Count','fontsize',14)
    title('Response Time','fontsize',14)

%cd('F:\Sam\DFD 2014 Presentation')
%print(gcf,'responseTime','-dpng','-r300')

t_c_data = xlsread('T_c.xlsx') ; 


pitchUpInd = find(t_c_data(:,3) >= 0) ; 
pitchDownInd = find(t_c_data(:,3) < 0) ; 

hcorrection = figure('paperPositionMode','auto') ;
    ylim = [18,52] ;
    xlim = [0 40] ;
    set(gca,'fontsize',14) ;
    hold on
    
    %q = quantile(abs(t_c_data(:,4)),[.25 .5 .75]) ; 
    meanTc = mean(t_c_data(:,4)) ;
    stdTc = std(t_c_data(:,4)) ;
    hline1 = plot([0 40], [meanTc-stdTc meanTc-stdTc], 'k--','LineWidth',1.5) ;
    hline2 = plot([0 40], [meanTc meanTc], 'k-','LineWidth',2.5) ;
    hline3 = plot([0 40], [meanTc+stdTc meanTc+stdTc], 'k--','LineWidth',1.5) ;
    
    hdata = plot(t_c_data(pitchUpInd,3), t_c_data(pitchUpInd,4), 'ko',...
        'MarkerSize', 8, 'MarkerFaceColor',[0 .7 0]) ; %,'MarkerEdgeColor',[0 .7 0]);
    plot(-t_c_data(pitchDownInd,3), t_c_data(pitchDownInd,4), 'ko',...
        'MarkerSize', 8, 'MarkerFaceColor',[0 .7 0]) ; %,'MarkerEdgeColor',[0 .7 0])
    
    %legend([hdata, hline],{'Data','Average'}, 'location','northwest')
    xlabel('\Delta Pitch [deg]','fontsize',14)
    ylabel('Correction Time, T_c [ms]','fontsize',14)
    title('Correction Time','fontsize',14)
    set(gca,'xlim',xlim)
    set(gca,'ylim',ylim)
    %print(gcf,'correctionTime','-dpng','-r300')
    
end    
    
    