%--------------------------------------------------------------------------
% function to fetch relevant data for PI controller fitting of pitch
% perturbations
%--------------------------------------------------------------------------
function [deltaPhiFront, c_pitch,fwdFlipTimes, backFlipTimes] = ...
    getPitchControllerData(data, debugFlag)

defineConstantsScript
smoothingParams = setSmoothingParams ; 

if isfield(data,'manualCorrRangeMS')
    manualCorrRangeMS_start = max([data.manualCorrRangeMS(1) -10]) ;
    manualCorrRangeMS_end = min([data.manualCorrRangeMS(2) 40]) ;
    manualCorrRangeMS = [manualCorrRangeMS_start, manualCorrRangeMS_end ] ;
else
    manualCorrRangeMS = [-10 40] ;
end
manualCorrRange = manualCorrRangeMS / 1000 ; 

%--------------------------------------------
% read in data
fwdFlipTimesR = data.fwdFlipTimesR ;
fwdFlipTimesL = data.fwdFlipTimesL ;
backFlipTimesR = data.backFlipTimesR ;
backFlipTimesL = data.backFlipTimesL ;

t = (data.params.startTrackingTime : data.params.endTrackingTime ) / 8000 ;

phiR = -data.anglesBodyFrame(:,PHIR) ;
phiL = data.anglesBodyFrame(:,PHIL) ;
bodyPitch = data.anglesLabFrame(:,BETA) ;

%--------------------------------------------
% restrict time range
correctedIndR = find(fwdFlipTimesR > manualCorrRange(1) & ...
    fwdFlipTimesR < manualCorrRange(2)) ;
correctedIndL = find(fwdFlipTimesL > manualCorrRange(1) & ...
    fwdFlipTimesL < manualCorrRange(2)) ;
correctedIndBackR = find(backFlipTimesR > manualCorrRange(1) & ...
    backFlipTimesR < manualCorrRange(2)) ;
correctedIndBackL = find(backFlipTimesL > manualCorrRange(1) & ...
    backFlipTimesL < manualCorrRange(2)) ;

fwdFlipTimesR = fwdFlipTimesR(correctedIndR) ;
fwdFlipTimesL = fwdFlipTimesL(correctedIndL) ;
backFlipTimesR = backFlipTimesR(correctedIndBackR) ;
backFlipTimesL = backFlipTimesL(correctedIndBackL) ;

%------------------------------------------
%spline smooth for stroke angle
phiEstErr = 1 ;
[sp_phiR, ~, ~] = mySplineSmooth(t(~isnan(phiR)),phiR(~isnan(phiR)),phiEstErr) ;
[sp_phiL, ~, ~] = mySplineSmooth(t(~isnan(phiL)),phiL(~isnan(phiL)),phiEstErr) ;

%------------------------------------------
%low-pass butterworth filter for body pitch data
pitch_filt = filterEulerAngle(bodyPitch,smoothingParams.pitch_filt_lvl) ;
c_pitch = fit(t',pitch_filt,'cubicinterp');
%pitch_smoothed = c_pitch(t) ;
%pitchVel_smoothed = differentiate(c_pitch, t) ; 
%c_pitchVel = fit(t',pitchVel_smoothed,'linearinterp');

%------------------------------------------
% if only using one wing
if ~isfield(data,'oneWing')
    phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
    phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;
elseif strcmp(data.oneWing,'L')
    phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;
    phiFrontR = phiFrontL ;
    fwdFlipTimesR = fwdFlipTimesL ; 
elseif strcmp(data.oneWing,'R')
    phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
    phiFrontL = phiFrontR ;
    fwdFlipTimesL = fwdFlipTimesR ; 
else 
    phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
    phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;
end
    
%--------------------------------------------------------------------------
% correct for array length differences
if length(fwdFlipTimesR) == length(fwdFlipTimesL)
    phiFront = (phiFrontR + phiFrontL) / 2 ;
    %phiFront_errorBar = (phiFrontR - phiFrontL) / 2 ;
    fwdFlipTimes = (fwdFlipTimesR + fwdFlipTimesL ) /2 ;
elseif length(fwdFlipTimesR) < length(fwdFlipTimesL)
    idx = zeros(length(fwdFlipTimesR),1) ;
    for q = 1:length(fwdFlipTimesR)
        [~,minInd] = min(abs(fwdFlipTimesL - fwdFlipTimesR(q))) ; 
        idx(q) = minInd ;
    end
    fwdFlipTimes = (fwdFlipTimesR + fwdFlipTimesL(idx)) / 2 ;
    phiFront = (phiFrontR + phiFrontL(idx)) / 2 ;
    %phiFront_errorBar = (phiFrontR + phiFrontL(idx)) / 2 ;
elseif length(fwdFlipTimesL) < length(fwdFlipTimesR)
    idx = zeros(length(fwdFlipTimesL),1) ;
    for q = 1:length(fwdFlipTimesL)
        [~,minInd] = min(abs(fwdFlipTimesR - fwdFlipTimesL(q))) ; 
        idx(q) = minInd ;
    end
    fwdFlipTimes = (fwdFlipTimesL + fwdFlipTimesR(idx)) / 2 ;
    phiFront = (phiFrontL + phiFrontR(idx)) / 2 ;
    %phiFront_errorBar = (phiFrontL - phiFrontR(idx)) / 2 ;
end

if length(backFlipTimesR) == length(backFlipTimesL)
    backFlipTimes = (backFlipTimesR + backFlipTimesL ) /2 ;
elseif length(backFlipTimesR) < length(backFlipTimesL)
    idx = zeros(length(backFlipTimesR),1) ;
    for q = 1:length(backFlipTimesR)
        [~,minInd] = min(abs(backFlipTimesL - backFlipTimesR(q))) ;
        idx(q) = minInd ;
    end
    backFlipTimes = (backFlipTimesR + backFlipTimesL(idx)) / 2 ;
    
elseif length(backFlipTimesL) < length(backFlipTimesR)
    idx = zeros(length(backFlipTimesL),1) ;
    for q = 1:length(backFlipTimesL)
        [~,minInd] = min(abs(backFlipTimesR - backFlipTimesL(q))) ;
        idx(q) = minInd ;
    end
    backFlipTimes = (backFlipTimesL + backFlipTimesR(idx)) / 2 ;
    %phiBack = (phiBackL + phiBackR(idx)) / 2 ;
end

%--------------------------------------------------------------------------
% make plot to test data?
if debugFlag
    t_min = min([backFlipTimes(1) fwdFlipTimes(1)]) ;
    t_max = max([backFlipTimes(end) fwdFlipTimes(end)]) ;
    t_range = linspace(t_min, t_max, 300) ; 
    
    pitch_smoothed = c_pitch(t_range) ;
    [pitchVel, pitchAccel] = differentiate(c_pitch, t_range) ;
    
    figure ;
    subplot(2,2,1)
    hold on
    plot(1000*t_range, pitch_smoothed,'b-','LineWidth',1.5)
    plot(1000*t, bodyPitch, 'k.')
    %plot(1000*t_peak, pitch_smoothed(ind_peak), 'ro', 'MarkerFaceColor','r')
    ylabel('\theta [deg]')
    
    subplot(2,2,3)
    plot(1000*t_range, pitchVel,'m-','LineWidth',1.5)
    ylabel('Pitch Vel [deg/s]')
    
    subplot(2,2,2)
    plot(1000*t_range, pitchAccel,'r-','LineWidth',1.5)
    hold on
    %plot(1000*t_peak, pitchAccel(ind_peak),'bx')
    %set(gca,'xlim',1000*[t_min t_max])
    ylabel('Pitch Accel [deg/s^2]')
    xlabel('Time [ms]')
    
    subplot(2,2,4)
    hold on
    plot(1000*t_range, fnval(sp_phiR,t_range), 'r-')
    plot(1000*t_range, fnval(sp_phiL,t_range), 'b-')
    plot(1000*fwdFlipTimes, phiFront,'kv','MarkerFaceColor', 'k')
    %plot(1000*backFlipTimes, phiBack,'k^','MarkerFaceColor', 'k')
    set(gca,'xlim', 1000*[t_min t_max])
    ylabel('\phi [deg]')
end

prePertFlipInd = ( fwdFlipTimes < 0 ) ;
phiFrontAvg = mean(phiFront(prePertFlipInd)) ;
deltaPhiFront = phiFront - phiFrontAvg ; 


end