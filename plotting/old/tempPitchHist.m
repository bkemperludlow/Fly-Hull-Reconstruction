cd('D:\Janelia Flies\Analysis\No Perturbation')
load('analyzedMovies')

ExprNum =   analyzedMovies(:,1)' ; 
MovNum =    analyzedMovies(:,2)';

runSize = length(ExprNum);
numLags = 140 ; %for autocorrelation

dataPath = 'D:\Janelia Flies\Analysis\No Perturbation\' ;

mb258c_imptnta = [];
mb258c_tnte = [] ;
mb204b_imptnta = [];
mb204b_tnte = [] ;

mb258c_imptntaMedian = [] ;
mb258c_tnteMedian = [] ;
mb204b_imptntaMedian = [];
mb204b_tnteMedian = [] ;

mb258c_imptntaStd = [] ;
mb258c_tnteStd = [] ;
mb204b_imptntaStd = [];
mb204b_tnteStd = [] ;

mb258c_imptntaACF = [] ;
mb258c_tnteACF = [] ;
mb204b_imptntaACF = [];
mb204b_tnteACF = [] ;

for i = 1:runSize
    if MovNum(i) < 10
        zstr = '00' ;
    elseif MovNum(i) < 100 ;
        zstr = '0' ;
    else
        zstr = '' ;
    end
    movNumStr = [zstr num2str(MovNum(i))] ;
    cd([dataPath 'Expr_' num2str(ExprNum(i)) '_mov_' movNumStr ])
    
    %load(['Expr' num2str(ExprNum(i)) 'mov' movNumStr '_results.mat']) ;
    try 
        load('axisHats.mat') ;
    catch
        try
            load(['Expr' num2str(ExprNum(i)) 'mov' movNumStr '_results.mat']) ;
            axisHats = data.AHat ;
        catch
            load(['Expr_' num2str(ExprNum(i)) '_mov_' movNumStr '_results.mat']) ;
            axisHats = data.AHat ;
        end
    end
    
    %bodyAxis = data.AHat ;
    %bodyPitch = (180/pi)*asin(bodyAxis(:,3)) ;
    bodyPitch = (180/pi)*asin(axisHats(:,3)) ;
    bodyPitchAvg = zeros(floor(length(bodyPitch)/8),1) ; 
    for j = 1:length(bodyPitchAvg)
        i1 = 8*(j-1)+1 ;
        i2 = 8*j ;
        bodyPitchAvg(j) = mean(bodyPitch(i1:i2)) ;
    end
    acf = nan(numLags+1,1) ;
    acfTemp = autocorr(bodyPitch,min([numLags, length(bodyPitch)-1])) ;
    acf(1:length(acfTemp)) = acfTemp(1:end) ;
    
    %added by SW temporarily
    %{
    indNeg = find(bodyPitch < -20);
    if ~isempty(indNeg) 
        disp('Experiment Number')
        disp(ExprNum(i))
        disp('Movie Number')
        disp(MovNum(i))
    end    
    %}
    
    if (ExprNum(i) == 2) || (ExprNum(i) == 6) || (ExprNum(i) == 8) || (ExprNum(i) == 9)
        mb258c_tnte = [mb258c_tnte; bodyPitchAvg] ;
        mb258c_tnteMedian = [ mb258c_tnteMedian ; median(bodyPitch) ] ;
        mb258c_tnteStd = [ mb258c_tnteStd ; std(bodyPitch) ] ;
        mb258c_tnteACF = [ mb258c_tnteACF , acf ] ;
    elseif (ExprNum(i) == 1) || (ExprNum(i) == 3) || (ExprNum(i) == 7) || (ExprNum(i) == 10)
        mb258c_imptnta = [mb258c_imptnta; bodyPitchAvg] ;
        mb258c_imptntaMedian = [ mb258c_imptntaMedian ; median(bodyPitch) ] ;
        mb258c_imptntaStd = [ mb258c_imptntaStd ; std(bodyPitch) ] ;
        mb258c_imptntaACF = [ mb258c_imptntaACF , acf ] ;
    elseif (ExprNum(i) == 4) || (ExprNum(i) == 12) || (ExprNum(i) == 13)
        mb204b_tnte = [mb204b_tnte; bodyPitchAvg] ;
        mb204b_tnteMedian = [ mb204b_tnteMedian ; median(bodyPitch) ] ;
        mb204b_tnteStd = [ mb204b_tnteStd ; std(bodyPitch) ] ;
        mb204b_tnteACF = [  mb204b_tnteACF , acf ] ;
    elseif (ExprNum(i) == 5) || (ExprNum(i) == 11)
        mb204b_imptnta = [mb204b_imptnta; bodyPitchAvg] ;
        mb204b_imptntaMedian = [ mb204b_imptntaMedian ; median(bodyPitch) ] ;
        mb204b_imptntaStd = [ mb204b_imptntaStd ; std(bodyPitch) ] ;
        mb204b_imptntaACF = [ mb204b_imptntaACF , acf ] ;
    else
        disp('error')
    end
    clear('axisHats')
end

%xbins = 10:80 ;
%set(gca,'fontsize',14) ;
%hold on
nbins = 30 ;
[f1, x1] = hist(mb258c_imptnta,nbins) ;
[f2, x2] = hist(mb258c_tnte,nbins) ;
[f3, x3] = hist(mb204b_imptnta,nbins) ;
[f4, x4] = hist(mb204b_tnte,nbins) ;
%{    
hresponse1 = figure('paperPositionMode','auto') ;
    %b = bar([x1', x2', x3', x4'], [f1'/trapz(x1,f1), f2'/trapz(x2,f2), f3'/trapz(x3,f3), f4'/trapz(x4,f4)],...
    %    'grouped') ;
    b1 = bar(x1, f1/trapz(x1,f1),'edgecolor','none','facecolor','b') ;
    hold on
    b2 = bar(x2,f2/trapz(x2,f2),'edgecolor','none','facecolor','r') ;
    %childHandle = get(b2,'Children');
    %set(childHandle,'FaceAlpha',0.7); 
    
    legend({'MB258C UAS-IMPTNT-A','MB258C UAS-TNT-E'}, 'location','northwest')
    
    %{
    h = findobj(gca,'Type','patch');
    set(h(1),'FaceColor',[0 0 .8]);
    set(h(1),'EdgeColor','w');
    set(h(2),'FaceColor',[.8 0 0]);
    set(h(2),'EdgeColor','w');
    
    set(gca, 'Xtick', 0:3) ;
    set(gca, 'xlim',[0 3]) ; 
    %}
    xlabel('Body Pitch Angle [deg]')
    ylabel('Probability')
    title('PDF of Body Pitch Angle')

hresponse2 = figure('paperPositionMode','auto') ;
    %b = bar([x1', x2', x3', x4'], [f1'/trapz(x1,f1), f2'/trapz(x2,f2), f3'/trapz(x3,f3), f4'/trapz(x4,f4)],...
    %    'grouped') ;
    b3 = bar(x3, f3/trapz(x3,f3),'edgecolor','none','facecolor','c') ;
    hold on
    b4 = bar(x4,f4/trapz(x4,f4),'edgecolor','none','facecolor','g') ;
    %childHandle = get(b2,'Children');
    %set(childHandle,'FaceAlpha',0.7); 
    
    legend({'MB204B UAS-IMPTNT-A','MB204B UAS-TNT-E'}, 'location','northwest')
    
    %{
    h = findobj(gca,'Type','patch');
    set(h(1),'FaceColor',[0 0 .8]);
    set(h(1),'EdgeColor','w');
    set(h(2),'FaceColor',[.8 0 0]);
    set(h(2),'EdgeColor','w');
    
    set(gca, 'Xtick', 0:3) ;
    set(gca, 'xlim',[0 3]) ; 
    %}
    xlabel('Body Pitch Angle [deg]')
    ylabel('Probability')
    title('PDF of Body Pitch Angle')
%}    
hbox = figure ;
    maxLength = max([length(mb258c_imptnta), length(mb258c_tnte), length(mb204b_imptnta),...
        length(mb204b_tnte)]) ;
    
    temp1 = padarray(mb258c_imptnta,[maxLength - length(mb258c_imptnta),0],nan,'post') ;
    temp2 = padarray(mb258c_tnte,[maxLength - length(mb258c_tnte),0],nan,'post') ;
    temp3 = padarray(mb204b_imptnta,[maxLength - length(mb204b_imptnta),0],nan,'post') ;
    temp4 = padarray(mb204b_tnte,[maxLength - length(mb204b_tnte),0],nan,'post') ;

    X = [temp1 , temp2, temp3, temp4] ;
    G = {'MB258C-IMPTNT-A', 'MB258C-TNT-E', 'MB204B-IMPTNT-A', 'MB204B-TNT-E'} ;
    
    allcolorgroups = {'A','B','A','B'};
    cmap = hsv2rgb([0.3 0.6 0.6; 0 0.6 0.6]);
    
    hold on
    boxplot(X,G,'plotstyle','compact','symbol','ko',...
        'colorgroup',allcolorgroups,'colors',cmap,'medianstyle','line',...
        'jitter',.1);
    % make median lines black and big
    %set(findobj(gcf,'Tag','Median'),'Color',[0 0 0],'LineWidth',2);

    % make outlier dots gray and big
    set(findobj(gcf,'Tag','Outliers'),'MarkerEdgeColor',[0.4 0.4 0.4],'MarkerSize',3);
    %boxplot(X,G,'plotstyle','compact','boxstyle','filled')

    
hcorr = figure ;
    hold on
    for i = 1:size(mb258c_imptntaACF,2)
        plot(mb258c_imptntaACF(:,i),'Color',[0 .7 0],'LineWidth',2)
    end
    for i = 1:size(mb258c_tnteACF,2)
        plot(mb258c_tnteACF(:,i),'Color',[.7 0 0],'LineWidth',2)
    end
    for i = 1:size(mb204b_imptntaACF,2)
        plot(mb204b_imptntaACF(:,i),'Color',[0 .7 0],'LineWidth',2)
    end
    for i = 1:size(mb204b_tnteACF,2)
        plot(mb204b_tnteACF(:,i),'Color',[.7 0 0],'LineWidth',2)
    end

hhist = figure('Position',[140 100 500 700])   ; 
    xlim = [-5 70] ;
    ylim = [0 0.07] ;
    
    s1 = subplot(4,1,1) ;
    b1 = bar(x1, f1,'edgecolor','k','facecolor',[0 .7 0]) ;
    %set(gca,'xlim',xlim)
    %set(gca,'ylim',ylim)
    set(gca,'fontsize',12)
    
    s2 = subplot(4,1,2) ;
    b2 = bar(x2, f2,'edgecolor','k','facecolor',[.7 0 0]) ;
    %set(gca,'xlim',xlim)
    %set(gca,'ylim',ylim)
    set(gca,'fontsize',12)
    
    s3 = subplot(4,1,3) ;
    b3 = bar(x3, f3,'edgecolor','k','facecolor',[0 .7 0]) ;
    %set(gca,'xlim',xlim)
    %set(gca,'ylim',ylim)
    set(gca,'fontsize',12)
    
    s4 = subplot(4,1,4) ;
    b4 = bar(x4, f4,'edgecolor','k','facecolor',[.7 0 0]) ;
    %set(gca,'xlim',xlim)
    %set(gca,'ylim',ylim)
    set(gca,'fontsize',12)
    
    
