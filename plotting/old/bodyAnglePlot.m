function h = bodyAnglePlot(data, pitchFlag, yawFlag, rollFlag, manyLinesFlag)

defineConstantsScript

yawColor1 = [153 85 255]/255  ;%[0 139 139]/255 ;
yawColor2 = [198 175 233]/255 ;
rollColor1 = [55 155 0]/255 ; %[176 48 176]/255 ;
rollColor2 = [175 233 175]/255 ;%[176 48 176]/255 + [-.2 .1 -.2]  ;
pitchColor1 = [200 55 55]/255 ; %[176 48 48]/255 ;%[1 .2 .2] ;
pitchColor2 = [222 135 135]/255 ; 
pertColor = [255 255 51]/255 ;

markersize = 4 ;
linewidth1 = 3 ;
linewidth2 = .5 ;

patchColor = [1 1 1 ] * 0.8 ; 
fwdFlipTimesR = data.fwdFlipTimesR ;
backFlipTimesR = data.backFlipTimesR ;

%% Load relevant things from data
if isfield(data,'t')
    t = data.t ;
else
    t1_temp = data.params.startTrackingTime/8000 ;
    t2_temp = data.params.endTrackingTime/8000 ;
    dt_temp = 1/data.params.fps ;
    t = t1_temp:dt_temp:t2_temp ;
end
tms = t*1000 ;

%xlim = [-10, tms(end)] ;
%ylim = [-22 28] ;

%rhoTimes = data.rhoTimes ;

bodyyaw_raw = data.anglesLabFrame(:,PSI) ;
YawEstErr = 1.2 ;
[~,bodyyaw,~] = mySplineSmooth(t,bodyyaw_raw,YawEstErr) ;

bodypitch_raw = data.anglesLabFrame(:,BETA) ;
PitchEstErr = .6 ;
[~,bodypitch,~] = mySplineSmooth(t,bodypitch_raw,PitchEstErr) ;

bodyroll = data.anglesLabFrame(:,RHO) ;

%% Define reference angle value
if isfield(data.params,'pulseLengthMS')
    %pertFlag = true ;
    
    pulseLength = data.params.pulseLengthMS ;
    %tsfvec = [0 pulseLength pulseLength 0 0] ;
    %patchColor = [1 1 1 ] * 0.8 ; 
    %faceAlpha = 1 ;

    zeroInd = find(t == 0) ;
    %pitchRef = bodypitch(zeroInd) ;
    %yawRef = bodyyaw(zeroInd) ;
    %rollRef = bodyroll(zeroInd) ;
else
    zeroInd = find(t == 0) ;
    %pertFlag = false ;
    %pitchRef = median(bodypitch) ;
    %yawRef = median(bodyyaw) ;
    %rollRef = median(bodyroll) ;
end

%% Plot many data lines, if necessary
if manyLinesFlag
    cd('F:\luca\Analysis\metadata\Delta Body Kinematics')
    load('pitchUpArray.mat')
    load('pitchDownArray.mat')
end
%% Make main figure 
h = figure ;
    hold on
    
    if rollFlag
        if manyLinesFlag
           for i = 1:size(pitchUpArray,1)
               plot(1000*pitchUpArray{i,3},pitchUpArray{i,4},'-','LineWidth',linewidth2,'Color',rollColor2)
           end
           for i = 1:size(pitchDownArray,1)
               plot(1000*pitchDownArray{i,3},pitchDownArray{i,4},'-','LineWidth',linewidth2,'Color',rollColor2)
           end
        end  
        %hroll = plot(tms(1:4:end),bodyroll(1:4:end),'ko','MarkerFaceColor',rollColor1,...
        %    'MarkerEdgeColor',rollColor1,'MarkerSize',markersize) ;
        plot(tms,bodyroll,'Color',rollColor1,'LineWidth',linewidth1) ;
    end
    if yawFlag
        if manyLinesFlag
           for i = 1:size(pitchUpArray,1)
               plot(1000*pitchUpArray{i,3},pitchUpArray{i,5},'-','LineWidth',linewidth2,'Color',yawColor2)
           end
           for i = 1:size(pitchDownArray,1)
               plot(1000*pitchDownArray{i,3},pitchDownArray{i,5},'-','LineWidth',linewidth2,'Color',yawColor2)
           end
        end  
        %hyaw = plot(tms,bodyyaw,'ko','MarkerFaceColor',yawColor1,...
        %    'MarkerEdgeColor',yawColor1,'MarkerSize',markersize) ;
       plot(tms,bodyyaw- mean(bodyyaw),'Color',yawColor1,'LineWidth',linewidth1) ;
    end
    if pitchFlag
        if manyLinesFlag
           for i = 1:size(pitchUpArray,1)
               plot(1000*pitchUpArray{i,3},pitchUpArray{i,6},'-','LineWidth',linewidth2,'Color',pitchColor2)
           end
           for i = 1:size(pitchDownArray,1)
               plot(1000*pitchDownArray{i,3},pitchDownArray{i,6},'-','LineWidth',linewidth2,'Color',pitchColor2)
           end
        end
        %hpitch = plot(tms,bodypitch,'ko','MarkerFaceColor',pitchColor1,...
        %    'MarkerEdgeColor',pitchColor1,'MarkerSize',markersize) ; 
        plot(tms,bodypitch,'Color',pitchColor1,'LineWidth',linewidth1) ;
    end
    %{
    if pertFlag
        ylim = get(gca,'ylim') ;
        avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
        hf = fill(tsfvec , avec,'y') ;
        uistack(hf,'bottom')
        faceAlpha = .5 ;
        set(hf,'facecolor',pertColor,'facealpha',faceAlpha) ;
        set(hf,'edgecolor','k') ;
        set(hf,'HandleVisibility','off')
    end
    %}
    %plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
    %set(gca,'xlim',xlim)
    %set(gca,'ylim',ylim)
    %xlabel('Time [ms]','fontsize',12)
    %ylabel('Angle [deg]','fontsize',12)
    %legend({'Roll','Yaw','Pitch'},'location','northwest')
    %title('Deviation in Body Angles from Median Values','fontsize',12)
end