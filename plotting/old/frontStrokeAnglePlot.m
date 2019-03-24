function [hfront, hback] = frontStrokeAnglePlot(data,legendFlag,ylim1,ylim2)
%frontStrokeAnglePlot(data,false,[15 50],[150 185])
%-----------------------------------------------------------------
defineConstantsScript

linewidth = 1 ;
ColorFront = [0.91 0.41 0.17] ;
ColorBack = [0.5859 0.4805 0.7109] ;
%lavender = [0.5859 0.4805 0.7109]
%deep carrot orange = [0.91 0.41 0.17]
patchColor = [1 1 1 ] * 0.8 ; 
faceAlpha = 1 ;

if isfield(data,'t')
    t = data.t ;
else
    t1_temp = data.params.startTrackingTime/8000 ;
    t2_temp = data.params.endTrackingTime/8000 ;
    dt_temp = 1/data.params.fps ;
    t = t1_temp:dt_temp:t2_temp ;
end
tms = t*1000 ;

if isfield(data.params,'pulseLengthMS')
    pertFlag = true ;
    pulseLength = data.params.pulseLengthMS ;
    tsfvec = [0 pulseLength pulseLength 0 0] ;
else
    pertFlag = false ;
end

if (isfield(data,'ignoreFrames'))
    ignoreFrames = data.ignoreFrames ;
else
    ignoreFrames = [] ;
end

if isfield(data, 'manualCorrRangeMS')
    xlim = data.manualCorrRangeMS ;
elseif isfield(data, 'correctionTime')
    xlim = data.correctionTime ;
else
    xlim = [tms(1) tms(end)] ;
end

fwdFlipTimesR = data.fwdFlipTimesR ;
backFlipTimesR = data.backFlipTimesR ;
%fwdFlipTimesL = data.fwdFlipTimesL ;
%backFlipTimesL = data.backFlipTimesL ;

phiR = -data.anglesBodyFrame(:,PHIR) ;
phiL =  data.anglesBodyFrame(:,PHIL) ;

fwdFlipInd = zeros(size(fwdFlipTimesR)) ;
backFlipInd = zeros(size(backFlipTimesR)) ;

for i = 1:length(fwdFlipInd)
    fwdFlipInd(i) = find(t == fwdFlipTimesR(i)) ;
end
for j = 1:length(backFlipInd)
    backFlipInd(j) = find(t == backFlipTimesR(j)) ;
end

phiFront = .5*(phiR(fwdFlipInd) + phiL(fwdFlipInd)) ;
phiBack = .5*(phiR(backFlipInd) + phiL(backFlipInd)) ;

%% Front stroke
hfront = figure ;

%set(gcf, 'PaperPositionMode', 'manual');
%set(gcf, 'PaperUnits', 'inches');
%set(gcf, 'PaperPosition', [2 1 9 3]);
hold on;
set(gca,'fontsize',10) ;

if pertFlag
    avec = [ylim1(1) ylim1(1) ylim1(2) ylim1(2) ylim1(1)] ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
end

%plot(fwdFlipTimesR*1000,phiFront,'k-') ;

plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
plot(fwdFlipTimesR*1000,phiFront,'-o','LineWidth',linewidth,...
    'Color', ColorFront,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor', ColorFront,...
    'MarkerSize',8) ;
%hback = plot(backFlipTimesR*1000,phiBack,'-ko','MarkerFaceColor', ColorBack,'LineWidth',linewidth) ;

set(gca,'ylim',ylim1)
set(gca,'xlim',xlim)

xlabel('Time [ms]','fontsize',10)
ylabel({'Front Stroke Angle'; '  [deg]'},'fontsize',10)
title('Front Stroke Angle vs Time','fontsize',10)

if legendFlag
   legend({'\phi^\text{w}_\text{front}'},'location','northwest')
   
end

%% Back stroke
hback = figure ;
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [2 1 9 3]);
hold on;
set(gca,'fontsize',10) ;

if pertFlag
    avec = [ylim2(1) ylim2(1) ylim2(2) ylim2(2) ylim2(1)] ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
end

%hfront = plot(fwdFlipTimesR*1000,phiFront,'-ko','MarkerFaceColor', ColorFront) ;%,'LineWidth',linewidth) ;
plot(backFlipTimesR*1000,phiBack,'-o','LineWidth',linewidth,...
    'Color', ColorBack,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor', ColorBack,...
    'MarkerSize',8) ;
plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
set(gca,'ylim',ylim2)
set(gca,'xlim',xlim)

xlabel('Time [ms]','fontsize',10)
ylabel({'Back Stroke Angle'; '  [deg]'},'fontsize',10)
title('Back Stroke Angle vs Time','fontsize',10)

if legendFlag
   legend({'\phi^\text{w}_\text{back}'},'location','northwest')
end

end
