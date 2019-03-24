%NEEDS TO BE CLEANED UP A TON. JUST USE MK2. ALSO NEEDS ADJUSTMENT FOR PITCH DOWN VS UP

csvfilename = 'pitchUpMetadata2.csv' ;

ExprNum = [ 7 7 7 7 7 7 7 7 7 7 7 7 7];
MovNum = {'008' '009' '015' '019' '023' '024' '029' ...
    '036' '045' '046' '054' '058' '062'};

if (length(ExprNum) ~= length(MovNum))
    disp('Fix either ExprNum or MovNum arrays; they do not match')
    return ;
end

runSize = length(ExprNum) ;

defineConstantsScript
patchColor = [1 1 1 ] * 0.8 ; 
faceAlpha = 1 ;

%want to ultimately make a table with everything, so need to initialize arrays
deltaPitch = zeros(runSize,1) ;
deltaVel = zeros(runSize,1) ;
deltaAccel = zeros(runSize,1) ;
PitchEstErr = .6*ones(runSize,1) ;
%windowSize = 30*ones(runSize,1) ;
Mov = zeros(runSize,1) ;
Expr = ExprNum' ; 


for i = 4:4
    
    %Load data
    datapath = strcat('F:\luca\Analysis\pitch up\Expr_',...
        num2str(ExprNum(i)), '_mov_',MovNum{i}) ;
    cd(datapath)
    
    datafilename = strcat(datapath,'\Expr',num2str(ExprNum(i)), ...
        'mov',MovNum{i}, '_Data_manually_corrected.mat') ;
    
    load(datafilename) ;
    Mov(i) = str2num(MovNum{i}) ;
    
    %Check ignoreFrames and correctionTime
    if (isfield(data,'ignoreFrames'))
        ignoreFrames = data.ignoreFrames ;
    else
        ignoreFrames = [] ;
    end
    
    if (isfield(data,'correctionTime'))
        correctionTime=data.correctionTime;
    elseif(isfield(data,'manualCorrRangeMS'))
        correctionTime=data.manualCorrRangeMS;
    end
    
    if isfield(data.params,'pulseLengthMS')
        pulseLength = data.params.pulseLengthMS ;
    elseif (~isfield(data.params,'pulseLengthMS')) && (ExprNum(i) == 7)
        pulseLength = 5.8 ; %ms
    else
        pulseLength = 8 ; %ms
    end
    
    %Get flip times
    flipTimesFlag = 0 ;
    if isfield(data,'fwdFlipTimesR') && isfield(data,'backFlipTimesR')
        fwdFlipTimesR = data.fwdFlipTimesR ; 
        backFlipTimesR = data.backFlipTimesR ; 
        flipTimesFlag = 1 ;
    end
    
    %Get the body pitch angle
    if (isfield(data, 'anglesBodyFrame'))
        bodyPitch = data.anglesBodyFrame(:,BETA) ;
    elseif (isfield(data,'AHat'))
        axisHats = data.AHat ;
        bodyPitch = (180/pi)*asin(axisHats(:,3)) ;
    else
        disp('Check data structure (before you wreck data structure)')
    end
    
    %exclude bad frames
    startTime = data.params.startTrackingTime ;
    endTime = data.params.endTrackingTime ;
    fps = data.params.fps ;
    t = (startTime:endTime)/fps ; %seconds
    tms = 1000*t ; %ms
    dt = 1/fps ;
    zeropoint = find(t == 0) ;
    
    bodyPitch(ignoreFrames) = NaN ; 

    ind = find(~isnan(bodyPitch)) ;
    currtvec = t(ind) ;
    currBodyPitch = bodyPitch(ind)' ;
    
  
    currBeforePertInd = find(currtvec <= 0) ;
    currPertInd = find(0 <= currtvec & currtvec <= pulseLength/1000) ;
    currAfterPertInd = find(currtvec >= pulseLength/1000) ;
    weights = ((currtvec(end)-currtvec(1))/length(currtvec))*ones(size(currtvec)) ;
    weights([1 end]) = .5*weights([1 end]) ;
    weights([currPertInd(1) currPertInd(end)]) = 500*weights([currPertInd(1) currPertInd(end)]) ; %Need to figure out what proper weights are
    
    %smooth the data
    
    [sp_bodyPitch1, ~, ~] = mySplineSmooth(currtvec(currBeforePertInd), currBodyPitch(currBeforePertInd), PitchEstErr(i), weights(currBeforePertInd)) ;
    [sp_bodyPitch2, ~, ~] = mySplineSmooth(currtvec(currPertInd), currBodyPitch(currPertInd), PitchEstErr(i), weights(currPertInd)) ;
    [sp_bodyPitch3, ~, ~] = mySplineSmooth(currtvec(currAfterPertInd), currBodyPitch(currAfterPertInd), PitchEstErr(i), weights(currAfterPertInd)) ;
    
    bodyPitch_sg = sgolayfilt(currBodyPitch,5,201) ;

    beforePert = find(t < 0) ;
    pert = find(0 <= t & t <= pulseLength/1000) ;
    afterPert = find(t > pulseLength/1000) ;
    
    bodyPitch_smooth = [fnval(sp_bodyPitch1,t(beforePert)), fnval(sp_bodyPitch2,t(pert)),...
        fnval(sp_bodyPitch3,t(afterPert))]';
    pitchVelocity = [fnval( fnder(sp_bodyPitch1, 1), t(beforePert)), ...
        fnval( fnder(sp_bodyPitch2, 1) ,t(pert)), ...
        fnval( fnder(sp_bodyPitch3, 1),t(afterPert))]';
    pitchAccel = [fnval( fnder(sp_bodyPitch1, 2), t(beforePert)), ...
        fnval( fnder(sp_bodyPitch2, 2) ,t(pert)), ...
        fnval( fnder(sp_bodyPitch3, 2),t(afterPert))]';
    
    %{
    [sp_pitch, ~ , ~] = mySplineSmooth(currtvec, currBodyPitch, PitchEstErr(i)) ;%, weights) ;
    bodyPitch_smooth = fnval(sp_pitch,t)' ;
    pitchVelocity = fnval( fnder(sp_pitch,1) ,t)' ;
    pitchAccel = fnval( fnder(sp_pitch,2) ,t)' ;
    %}
    %Extract some of the relevant quantities (max change in accel, max velocity, etc.
    
    maxPitchVel = max(pitchVelocity) ;
    minPitchVel = min(pitchVelocity) ;
    baselineVel = mean(pitchVelocity(1:zeropoint)) ;
    deltaVel(i) = maxPitchVel - baselineVel ; 
    
    maxPitchAccel = max(pitchAccel) ;
    minPitchAccel = min(pitchAccel) ;
    baselineAccel = mean(pitchAccel(1:zeropoint)) ;
    deltaAccel(i) = abs(minPitchAccel - baselineAccel) ; 
    
    %Define region for perturbation rectangle. 
    tsfvec = [0 pulseLength pulseLength 0 0] ;
    
    ylimPitch = [(minPitch-5) (maxPitch+5)] ;
    ylimVel = [(minPitchVel-abs(.2*minPitchVel)) (maxPitchVel+abs(.2*maxPitchVel))] ;
    ylimAccel = [(minPitchAccel-abs(.2*minPitchAccel)) (maxPitchAccel+abs(.2*minPitchAccel))] ;
    xlim = [tms(1) tms(end)] ;
    
    avec1 = [ylimPitch(1) ylimPitch(1) ylimPitch(2) ylimPitch(2) ylimPitch(1)] ;
    avec2 = [ylimVel(1) ylimVel(1) ylimVel(2) ylimVel(2) ylimVel(1)] ;
    avec3 = [ylimAccel(1) ylimAccel(1) ylimAccel(2) ylimAccel(2) ylimAccel(1)] ;
    
    %plot to check that splines make sense. Needs to have perturbation and wing- 
    %stroke background probably
    %{
    hpolyfit = figure ; 
        plot(1000*currtvec,currBodyPitch,'kx') 
        ylim = get(gca,'ylim') ;
        box on; grid on;
        hold on
        plot(1000*t,bodyPitch_poly(1)*t.^2 + bodyPitch_poly(2)*t +bodyPitch_poly(3)*ones(size(t)),'LineWidth',2)

        plot(1000*currtvec(polyInd(1)),currBodyPitch(polyInd(1)),'r.','MarkerSize',20)
        plot(1000*currtvec(polyInd(end)),currBodyPitch(polyInd(end)),'r.','MarkerSize',20)

        plot(1000*t,bodyPitch_smooth, 'g--','LineWidth',2)

        plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
        set(gca, 'ylim',ylim)
    %}
    hpitch = figure ;
    subplot(3,1,1)
        hf = fill(tsfvec , avec1,'y') ;
        set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
        hold on;
        plot(tms, bodyPitch_smooth, 'r-','LineWidth',2)
        %plot(currtvec*1000,bodyPitch_sg,'b--','LineWidth',1.5)
        plot(tms,bodyPitch,'kx','MarkerSize',4.5)
        box on ; grid on;
        ylabel('Body Pitch Angle [deg]')
        title(strcat('Expr ',num2str(ExprNum(i)),' Mov ',MovNum(i)))
        if flipTimesFlag
            plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
        end
        set(gca,'ylim',ylimPitch) ;
        set(gca,'xlim',xlim) ;
        
    subplot(3,1,2)
        hf = fill(tsfvec , avec2,'y') ;
        set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
        hold on;
        plot(tms, pitchVelocity, 'Color', [.5 0 .5], 'LineWidth',2)
        box on ; grid on;
        ylabel('Pitch Velocity [deg/s]')
        if flipTimesFlag
            plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
        end
        set(gca,'ylim',ylimVel) ;
        set(gca,'xlim',xlim) ;
        
    subplot(3,1,3)
        hf = fill(tsfvec , avec3,'y') ;
        set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
        hold on;
        plot(tms, pitchAccel, 'Color', [0 .5 .5], 'LineWidth',2)
        
        box on ; grid on;
        ylabel('Pitch Acceleration[deg/s^2]')
        xlabel('Time [ms]')
        if flipTimesFlag
            plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
        end
        set(gca,'ylim',ylimAccel) ;
        set(gca,'xlim',xlim) ;
    
end
%{
cd('F:\luca\Analysis\metadata') ;
currmetadata = table(Expr, Mov, deltaPitch, deltaVel, deltaAccel, ...
    PitchEstErr) ; 

if (exist(csvfilename) == 2)
    oldmetadata = readtable(csvfilename) ;
    metadata = union(oldmetadata,currmetadata) ;
    writetable(metadata, csvfilename) ; 
else
    metadata = currmetadata ; 
    writetable(metadata, csvfilename) ;
end
%}
