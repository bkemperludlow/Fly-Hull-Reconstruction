function h = wingAnglePlot(data,angleType,dataFlag,legendFlag,ylim)
%h = wingAnglePlot(data,0,true,true,[0 180])
%h = wingAnglePlot(data,1,true,true,[0 180])
%h = wingAnglePlot(data,2,true,true,ylim)
%--------------------------------------------------------------------------
%Produces a figure of some angle variables for a single movie
%
%Inputs:
%   data - typical data structure
%   angleType - integer that indicates which type of angle to plot
%       0 = stroke angle
%       1 = wing pitch
%       2 = deviation angle
%   dataFlag - plot real data or just spline
%   legendFlag - legend or no
%   ylim - limits of the y axis
%Output:
%   h - figure 
%--------------------------------------------------------------------------
defineConstantsScript

linewidth = .75 ;
ColorL = [.2 .2 1] ;
ColorR = [1 0 0] ;
markersize = 3;

bodyPitchFlag = false ;

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
%{
fwdFlipTimesRInd = zeros(size(fwdFlipTimesR)) ;
backFlipTimesRInd = zeros(size(backFlipTimesR)) ;

for j = 1:length(fwdFlipTimesR)
    [~,tmp] = min(abs(t - fwdFlipTimesR(j))) ;
    fwdFlipTimesRInd(j) = tmp ;
end
for k = 1:length(backFlipTimesR)
    [~,tmp] = min(abs(t - backFlipTimesR(k))) ;
    backFlipTimesRInd(k) = tmp ;
end
onlyBackStrokes = [] ;
for it = 1:length(fwdFlipTimesR)
    onlyBackStrokes = [onlyBackStrokes, fwdFlipTimesRInd(it):backFlipTimesRInd(it+1)] ;
end
onlyFwdStrokes = [] ;
for it2 = 1:length(fwdFlipTimesR)
    onlyFwdStrokes = [onlyFwdStrokes, backFlipTimesRInd(it2):fwdFlipTimesRInd(it2)] ;
end
%fwdFlipTimesL = data.fwdFlipTimesL ;
%backFlipTimesL = data.backFlipTimesL ;
%}

switch angleType
    case 0
        angleR = -data.anglesBodyFrame(:,PHIR) ;
        angleL = data.anglesBodyFrame(:,PHIL) ;
        EstErr = 1 ;
        labelStr = 'Wing Stroke Angle' ;
        
   case 1
        angleR = (180/pi)*mod(unwrap((pi/180)*data.anglesBodyFrame(:,ETAR)),2*pi) ;
        angleL = (180/pi)*mod(unwrap((pi/180)*data.anglesBodyFrame(:,ETAL)),2*pi) ;
        EstErr = 2.5 ;
        labelStr = 'Wing Pitch Angle' ;
        
    case 2
        angleR = data.anglesBodyFrame(:,THETAR) ;
        angleL = data.anglesBodyFrame(:,THETAL) ;
        EstErr = 1 ;
        labelStr = 'Wing Deviation Angle' ; 
end

if (~isempty(ignoreFrames))
    angleR(ignoreFrames) = NaN ;
    angleL(ignoreFrames) = NaN ;
end

indR = find(~isnan(angleR)) ;
indL = find(~isnan(angleL)) ;
currtvecR = t(indR) ;
currtvecL = t(indL) ;
currangleR = angleR(indR) ;
currangleL = angleL(indL) ;

[sp_angleR, ~, ~] = mySplineSmooth(currtvecR,currangleR,EstErr) ;
[sp_angleL, ~, ~] = mySplineSmooth(currtvecL,currangleL,EstErr) ;

h = figure ;
%set(gcf, 'PaperPositionMode', 'manual');
%set(gcf, 'PaperUnits', 'inches');
%set(gcf, 'PaperPosition', [2 1 9 3]);
hold on;
set(gca,'fontsize',10) ;

if pertFlag
    avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 238 170]/255,'facealpha',faceAlpha) ;
end

if dataFlag
    hdataL = plot(t*1000,angleL,'o','MarkerFaceColor','none','MarkerEdgeColor',ColorL,'MarkerSize',markersize) ;
    hdataR = plot(t*1000,angleR,'o','MarkerFaceColor','none','MarkerEdgeColor',ColorR,'MarkerSize',markersize) ;
end


hR_sp = plot(t*1000,fnval(sp_angleR,t),'Color', ColorR,'LineWidth',linewidth) ;
hL_sp = plot(t*1000,fnval(sp_angleL,t),'Color', ColorL,'LineWidth',linewidth) ;
set(gca,'ylim',ylim)
set(gca,'xlim',xlim)
if bodyPitchFlag
    bodyPitch = data.anglesLabFrame(:,BETA) ;
    [~, bodyPitch_smooth,~] = mySplineSmooth(t, bodyPitch, 1) ;
    hbodyPitch = plot(t*1000,bodyPitch_smooth,'Color', 'k','LineWidth',linewidth,...
        'LineStyle','--') ;
end
plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
%plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, [1 1 1], false);
xlabel('Time [ms]') ;
ylabel([labelStr ' [deg]']) ;
title([labelStr ' vs Time'],'fontsize',10)


if legendFlag
    try
        legend([hf hL_sp hR_sp], {'Pert.', 'Left','Right'},'location','northwest')
    catch
        legend([hL_sp hR_sp], {'Left','Right'},'location','northwest')
    end
end

end
        
        
        
        
        
        
