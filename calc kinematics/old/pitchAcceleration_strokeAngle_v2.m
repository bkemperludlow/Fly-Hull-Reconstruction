rootPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\' ;
%cd(rootPath)

logPath = [rootPath 'Pitch Controller Analysis'] ;
cd(logPath)
dataLog = importdata('Pitch_Controller_Log.xlsx') ;

ExprNumArray = dataLog.data.Finished(:,1) ;
MovNumArray = dataLog.data.Finished(:,2) ;
pertTypeArray = dataLog.data.Finished(:,3) ;

ExprNumArray = ExprNumArray(~isnan(ExprNumArray)) ; 
MovNumArray = MovNumArray(~isnan(MovNumArray)) ; 
pertTypeArray = pertTypeArray(~isnan(pertTypeArray)) ; 

defineConstantsScript

accelerationAndPhiStruct = struct('ExprNum', [], 'MovNum', [], 't',[],'bodyPitch',[],...
    'pitchVelocity',[],'pitchAcceleration',[],'fwdFlipTimes', [], 'backFlipTimes', [],...
    'phiFront',[],'phiBack',[],'phiFrontPrePert',[],'phiBackPrePert',[],'maxPitchAccel',[],...
    'maxPitchVel',[], 'maxDeltaPitch', [],'deltaPhiFrontPeak',[],'deltaPhiBackPeak',[],...
    'deltaPitch',[],'deltaPhiFront',[],'peakInd',[],'pitchType',[],'flyType',[]);

d1 = designfilt('lowpassiir','FilterOrder',8,'SampleRate',8000, ...
        'HalfPowerFrequency',100,'DesignMethod','butter'); %hpf = 100
    
%cc = 1 ; 

for i = 1:length(MovNumArray)
    % load data 
    if pertTypeArray(i) == 1
        cd([rootPath 'Pitch Up']) 
    elseif pertTypeArray(i) == -1
        cd([rootPath 'Pitch Down']) 
    else 
        keyboard ; 
    end 
    
    ExprNum = ExprNumArray(i) ;
    MovNum = MovNumArray(i) ;

    if MovNum < 10
        zstr = '00' ;
    elseif MovNum < 100 ;
        zstr = '0' ;
    else
        zstr = '' ;
    end

    folderName = ['Expr_' num2str(ExprNum) '_mov_' zstr num2str(MovNum) ] ;
    fitFilename = ['Expr' num2str(ExprNum) 'mov' num2str(MovNum) '_controllerFit'] ;
    dataFilename1 = ['Expr' num2str(ExprNum) 'mov' zstr num2str(MovNum) '_Data_manually_corrected.mat'] ;
    dataFilename2 = ['Expr' num2str(ExprNum) 'mov' zstr num2str(MovNum) '_Data_manually_corrected_onlyPhiFront.mat'] ;
    dataFilename3 = ['Expr_' num2str(ExprNum) '_mov_' zstr num2str(MovNum) '_test.mat'] ;
    
    cd(folderName) 
    load( fitFilename )
    if exist(dataFilename1,'file') == 2
        load(dataFilename1)
    elseif exist(dataFilename2,'file') == 2
        load(dataFilename2)
    elseif exist(dataFilename3,'file') == 2
        load(dataFilename3)
    else
        disp('No data file')
        continue ;
    end
    
    if strcmp({controller_fit_struct.flyType},'control')
        flyType = 2 ;
    elseif strcmp({controller_fit_struct.flyType},'experimental') 
        flyType = 1 ;
    else 
        flyType = nan ;
    end
    
    % get relevant variables
    %if isfield(data,'manualCorrRangeMS')
    %    manualCorrRangeMS = data.manualCorrRangeMS ;
    %else
    %    manualCorrRangeMS = [-10 30] ;
    %end
    manualCorrRange = [controller_fit_struct.t(1), controller_fit_struct.t(end)]  / 1000 ;
    
    fwdFlipTimesR = data.fwdFlipTimesR ;
    fwdFlipTimesL = data.fwdFlipTimesL ;
    backFlipTimesR = data.backFlipTimesR ;
    backFlipTimesL = data.backFlipTimesL ;
    
    correctedIndFwdR = find(fwdFlipTimesR > manualCorrRange(1) & fwdFlipTimesR < manualCorrRange(2)) ;
    correctedIndFwdL = find(fwdFlipTimesL > manualCorrRange(1) & fwdFlipTimesL < manualCorrRange(2)) ;
    correctedIndBackR = find(backFlipTimesR > manualCorrRange(1) & backFlipTimesR < manualCorrRange(2)) ;
    correctedIndBackL = find(backFlipTimesL > manualCorrRange(1) & backFlipTimesL < manualCorrRange(2)) ;
    
    fwdFlipTimesR = fwdFlipTimesR(correctedIndFwdR) ;
    fwdFlipTimesL = fwdFlipTimesL(correctedIndFwdL) ;
    backFlipTimesR = backFlipTimesR(correctedIndBackR) ;
    backFlipTimesL = backFlipTimesL(correctedIndBackL) ;
    
    if length(fwdFlipTimesR) == length(controller_fit_struct.fwdFlipTimes)
        if sum(abs(fwdFlipTimesR - controller_fit_struct.fwdFlipTimes)) < 0.00001
            data.oneWing = 'R' ; 
        end
    elseif length(fwdFlipTimesL) == length(controller_fit_struct.fwdFlipTimes)
        if sum(abs(fwdFlipTimesL - controller_fit_struct.fwdFlipTimes)) < 0.00001
            data.oneWing = 'L' ;
        end
    end
    
    t = (data.params.startTrackingTime : data.params.endTrackingTime ) / 8000 ; 

    phiR = -data.anglesBodyFrame(:,PHIR) ;
    phiL = data.anglesBodyFrame(:,PHIL) ;
    bodyPitch = data.anglesLabFrame(:,BETA) ;
    
    phiEstErr = 1 ;
    [sp_phiR, ~, ~] = mySplineSmooth(t,phiR,phiEstErr) ;
    [sp_phiL, ~, ~] = mySplineSmooth(t,phiL,phiEstErr) ;
    
    
    %low-pass butterworth filter for data
    
    pitch_filt = filtfilt(d1,bodyPitch) ;
    c_pitch = fit(t',pitch_filt,'cubicinterp');
    pitch_smoothed = c_pitch(t) ;
    [pitchVel, pitchAccel] = differentiate(c_pitch, t) ; 
    bodyPitch_prePert = mean(pitch_smoothed(t <= 0 & t > -0.01)) ; %pitch_temp(t0_ind) ;
    deltaPitch = abs(pitch_smoothed - bodyPitch_prePert) ; 
    [pks, locs] = findpeaks(deltaPitch(t > 0.007)) ;
    postPulse_ind = find(t  > 0.007, 1, 'first') ;
    ind_peak = locs(1) + postPulse_ind ;
    
    %check that the sign is right
    deltaPitch_signed = pitch_smoothed(ind_peak) - bodyPitch_prePert ; 
    if sign(deltaPitch_signed) ~= sign( pertTypeArray(i) ) 
        disp('Peak does not match pert type')
        keyboard ;
    end
    t_peak = t(ind_peak) ; 
    maxDeltaPitch = pitch_smoothed(ind_peak) - bodyPitch_prePert ; 
    maxPitchAccel = pitchAccel(ind_peak) ; 
    
    [~, pitchVelPeakLocs] = findpeaks(abs(pitchVel(t > 0))) ; 
    t0ind = find(t == 0) ; 
    maxPitchVel = pitchVel(pitchVelPeakLocs(1)+t0ind) ;
    
    %pitchEstErr = .48 ; %.5
    %[sp_pitch, ~,~] = mySplineSmooth(t,bodyPitch,pitchEstErr) ;
    %pitch_smoothed = fnval(sp_pitch, t) ;
    
    if ~isfield(data,'oneWing')
        phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
        phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;
        phiBackR = fnval(sp_phiR,backFlipTimesR) ;
        phiBackL = fnval(sp_phiL,backFlipTimesL) ;
    elseif strcmp(data.oneWing,'L')
        phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;
        phiFrontR = phiFrontL ;
        fwdFlipTimesR = fwdFlipTimesL ;
        phiBackL = fnval(sp_phiL,backFlipTimesL) ;
        phiBackR = phiBackL ;
        backFlipTimesR = backFlipTimesL ;
    elseif strcmp(data.oneWing,'R')
        phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
        phiFrontL = phiFrontR ;
        fwdFlipTimesL = fwdFlipTimesR ;
        phiBackR = fnval(sp_phiR,backFlipTimesR) ;
        phiBackL = phiBackR ;
        backFlipTimesL = backFlipTimesR ;
    else
        phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
        phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;
        phiBackR = fnval(sp_phiR,backFlipTimesR) ;
        phiBackL = fnval(sp_phiL,backFlipTimesL) ;
    end
    
    
    if length(fwdFlipTimesR) == length(fwdFlipTimesL)
        phiFront = (phiFrontR + phiFrontL) / 2 ;
        fwdFlipTimes = (fwdFlipTimesR + fwdFlipTimesL ) /2 ;
    elseif length(fwdFlipTimesR) < length(fwdFlipTimesL)
        idx = zeros(length(fwdFlipTimesR),1) ;
        for q = 1:length(fwdFlipTimesR)
            [~,minInd] = min(abs(fwdFlipTimesL - fwdFlipTimesR(q))) ;
            idx(q) = minInd ;
        end
        fwdFlipTimes = (fwdFlipTimesR + fwdFlipTimesL(idx)) / 2 ;
        phiFront = (phiFrontR + phiFrontL(idx)) / 2 ;
    elseif length(fwdFlipTimesL) < length(fwdFlipTimesR)
        idx = zeros(length(fwdFlipTimesL),1) ;
        for q = 1:length(fwdFlipTimesL)
            [~,minInd] = min(abs(fwdFlipTimesR - fwdFlipTimesL(q))) ;
            idx(q) = minInd ;
        end
        fwdFlipTimes = (fwdFlipTimesL + fwdFlipTimesR(idx)) / 2 ;
        phiFront = (phiFrontL + phiFrontR(idx)) / 2 ;
    end
    
    if length(backFlipTimesR) == length(backFlipTimesL)
        phiBack = (phiBackR + phiBackL) / 2 ;
        backFlipTimes = (backFlipTimesR + backFlipTimesL ) /2 ;
    elseif length(backFlipTimesR) < length(backFlipTimesL)
        idx = zeros(length(backFlipTimesR),1) ;
        for q = 1:length(backFlipTimesR)
            [~,minInd] = min(abs(backFlipTimesL - backFlipTimesR(q))) ;
            idx(q) = minInd ;
        end
        backFlipTimes = (backFlipTimesR + backFlipTimesL(idx)) / 2 ;
        phiBack = (phiBackR + phiBackL(idx)) / 2 ;
    elseif length(backFlipTimesL) < length(backFlipTimesR)
        idx = zeros(length(backFlipTimesL),1) ;
        for q = 1:length(backFlipTimesL)
            [~,minInd] = min(abs(backFlipTimesR - backFlipTimesL(q))) ;
            idx(q) = minInd ;
        end
        backFlipTimes = (backFlipTimesL + backFlipTimesR(idx)) / 2 ;
        phiBack = (phiBackL + phiBackR(idx)) / 2 ;
    end
    
    phiBackPrePert = mean(phiBack(backFlipTimes > -0.01 & backFlipTimes < 0)) ; 
    phiFrontPrePert = mean(phiFront(fwdFlipTimes > -0.01 & fwdFlipTimes < 0)) ;
    prePeakBackInd = find(backFlipTimes < t_peak,1,'last') ;
    prePeakFrontInd = find(fwdFlipTimes < t_peak,1,'last') ;
    
    deltaPhiBackPeak = phiBack(prePeakBackInd) - phiBackPrePert ; 
    deltaPhiFrontPeak = phiFront(prePeakFrontInd) - phiFrontPrePert ; 
    deltaPhiFront = phiFront - phiFrontPrePert ;
    deltaPhiBack = phiBack - phiBackPrePert ;
    
    if (0)
        t_min = 1000*min([backFlipTimes(1) fwdFlipTimes(1)]) ;
        t_max = 1000*max([backFlipTimes(end) fwdFlipTimes(end)]) ;
        
        figure ;
        subplot(2,2,1)
        hold on
        plot(1000*t, pitch_smoothed,'b-','LineWidth',1.5)
        plot(1000*t, bodyPitch, 'k.')
        plot(1000*t_peak, pitch_smoothed(ind_peak), 'ro', 'MarkerFaceColor','r')
        ylabel('\theta [deg]')
        
        subplot(2,2,3)
        plot(1000*t, pitchVel,'m-','LineWidth',1.5)
        ylabel('Pitch Vel [deg/s]')
        
        subplot(2,2,2)
        plot(1000*t, pitchAccel,'r-','LineWidth',1.5)
        hold on
        plot(1000*t_peak, pitchAccel(ind_peak),'bx')
        set(gca,'xlim',[t_min t_max])
        ylabel('Pitch Accel [deg/s^2]')
        xlabel('Time [ms]')
        
        subplot(2,2,4)
        hold on
        plot(1000*t, fnval(sp_phiR,t), 'r-')
        plot(1000*t, fnval(sp_phiL,t), 'b-')
        plot(1000*fwdFlipTimes, phiFront,'kv','MarkerFaceColor', 'k')
        plot(1000*backFlipTimes, phiBack,'k^','MarkerFaceColor', 'k')
        set(gca,'xlim',[t_min t_max])
        ylabel('\phi [deg]')
        %disp('blah') 
    end
    
    accelerationAndPhiStruct(i).ExprNum = ExprNum ;
    accelerationAndPhiStruct(i).MovNum = MovNum ;
    accelerationAndPhiStruct(i).t = t ;
    accelerationAndPhiStruct(i).bodyPitch = pitch_smoothed ;
    accelerationAndPhiStruct(i).pitchVelocity = pitchVel ;
    accelerationAndPhiStruct(i).pitchAcceleration = pitchAccel ;
    accelerationAndPhiStruct(i).fwdFlipTimes = fwdFlipTimes ;
    accelerationAndPhiStruct(i).backFlipTimes = backFlipTimes ;
    accelerationAndPhiStruct(i).phiFront = phiFront ;
    accelerationAndPhiStruct(i).phiBack = phiBack ;
    accelerationAndPhiStruct(i).phiFrontPrePert = phiFrontPrePert ;
    accelerationAndPhiStruct(i).phiBackPrePert = phiBackPrePert ;
    accelerationAndPhiStruct(i).deltaPhiFront = deltaPhiFront ;
    accelerationAndPhiStruct(i).deltaPhiBack = deltaPhiBack ;
    accelerationAndPhiStruct(i).maxPitchAccel = maxPitchAccel ;
    accelerationAndPhiStruct(i).maxPitchVel = maxPitchVel ;
    accelerationAndPhiStruct(i).maxDeltaPitch = maxDeltaPitch ;
    accelerationAndPhiStruct(i).deltaPhiFront = deltaPhiFrontPeak ;
    accelerationAndPhiStruct(i).deltaPhiBack = deltaPhiBackPeak ;
    accelerationAndPhiStruct(i).deltaPitch = deltaPitch ;
    accelerationAndPhiStruct(i).peakInd = ind_peak ;
    accelerationAndPhiStruct(i).flyType = flyType ;
    accelerationAndPhiStruct(i).pertType = pertTypeArray(i) ; 
    
    %cc = cc + 1 ;
    clear data 
end



%{
cd(logPath)
save accelerationAndPhiStruct accelerationAndPhiStruct
%}