function h = frontStrokeAnglePlot_v2(data,ylim1,ylim2)
%frontStrokeAnglePlot(data,[15 50],[150 185])
%-----------------------------------------------------------------
defineConstantsScript

linewidth = 3 ;
ColorFront = [0.91 0.41 0.17] ;
ColorBack = [0.5859 0.4805 0.7109] ;
%lavender = [0.5859 0.4805 0.7109]
%deep carrot orange = [0.91 0.41 0.17]
patchColor = [1 1 1 ] * 0.8 ; 
faceAlpha = 1 ;

%ylim1 = [15 50]; 
%ylim2 = [155 190] ;
%legendFlag = false ;

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


phiFront_smooth = csapi(fwdFlipTimesR,phiFront) ;
phiBack_smooth = csapi(backFlipTimesR,phiBack) ;
currtFront = fwdFlipTimesR(1):0.0001:fwdFlipTimesR(end) ;
currtBack = backFlipTimesR(1):0.0001:backFlipTimesR(end) ;

%% Make plot
h = figure ;
hold on;
if pertFlag
    avec = [ylim1(1) ylim1(1) ylim1(2) ylim1(2) ylim1(1)] ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
end

%plot(fwdFlipTimesR*1000,phiFront,'k-') ;

%plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
%[ax, hfront, hback] = plotyy(currtFront*1000,fnval(phiFront_smooth,currtFront),...
%    currtBack*1000,fnval(phiBack_smooth,currtBack)) ;
plot(currtFront*1000,fnval(phiFront_smooth,currtFront),'k-')
plot(fwdFlipTimesR*1000,phiFront,'ko','MarkerFaceColor',ColorFront,'MarkerSize',4) ;

%{
'-o','LineWidth',linewidth,...
    'Color', ColorFront,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor', ColorFront,...
    'MarkerSize',8) ;
%}
%set(gca,'xlim',xlim)
set(gca,'ylim',ylim1)
%{
set(hfront, 'LineWidth',linewidth) ;
set(hfront, 'Color',ColorFront) ;
set(ax(1), 'YColor', ColorFront)
set(ax(1),'ylim',ylim1)
set(ax(1), 'YTick',ylim1(1):10:ylim1(2))
set(ax(1),'xlim',xlim)
ylabel(ax(1), 'Front Stoke Angle [deg]')

set(hback, 'LineWidth',linewidth) ;
set(hback, 'Color',ColorBack) ;
set(ax(2), 'YColor', ColorBack)
set(ax(2),'ylim',ylim2)
set(ax(2), 'YTick',ylim2(1):10:ylim2(2))
ylabel(ax(2), 'Back Stoke Angle [deg]')
set(ax(2), 'YTickLabel',cellstr(strsplit(num2str(ylim2(1):10:ylim2(2)))))
set(ax(2),'xlim',xlim)
%}
%xlabel('Time [ms]','fontsize',10)
%ylabel({'Front Stroke Angle'; '  [deg]'},'fontsize',10)
%title('Front Stroke Angle vs Time','fontsize',10)

%if legendFlag
%   legend({'\phi^\text{w}_\text{front}'},'location','northwest')
%   
%end



end
