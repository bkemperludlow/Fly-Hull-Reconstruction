%--------------------------------------------------------------------------
% function to fetch relevant data for PI controller fitting of pitch
% perturbations
%--------------------------------------------------------------------------
function [deltaPhiFront, c_pitch, c_vel, fwdFlipTimes, backFlipTimes, steady_t] = ...
    getPitchControllerData_v3(data, debugFlag, steadyPhiFrontFlag)
% -----------------------
%% inputs and params
if ~exist('debugFlag','var')
    debugFlag = false ; 
end
if ~exist('steadyPhiFrontFlag','var')
    steadyPhiFrontFlag = true ; 
end

defineConstantsScript
smoothingParams = setSmoothingParams ; 

% other params:
noMatchCost = 1e-3 ; % (sec) used to solve matching problem for flip times
phiEstErr = 1 ; % (deg) estimated error for stroke angle, used for spline
fitSteadyFlag = true ; % if true, find steady-state phi front by linear fit. otherwise find sinle wingstroke with minimal pitch torque
debugFitFlag = true ; % make plots for steady-state estimate?
pulseDuration = data.pulseDuration ; % in milliseconds

if isfield(data,'manualCorrRangeMS')
    manualCorrRangeMS_start = max([data.manualCorrRangeMS(1) -20]) ;
    manualCorrRangeMS_end = min([data.manualCorrRangeMS(2) 40]) ;
    manualCorrRangeMS = [manualCorrRangeMS_start, manualCorrRangeMS_end ] ;
else
    manualCorrRangeMS = [-10, 40] ; %[-10, 40] %[-70 70] 
end
manualCorrRange = manualCorrRangeMS / 1000 ; 
%manualCorrRange = [-10,  60] / 1000 ; 
%--------------------------------------------
%% read in data
fwdFlipTimesR = data.fwdFlipTimesR ;
fwdFlipTimesL = data.fwdFlipTimesL ;
backFlipTimesR = data.backFlipTimesR ;
backFlipTimesL = data.backFlipTimesL ;

t = (data.params.startTrackingTime : data.params.endTrackingTime ) / 8000 ;
dt = nanmean(diff(t)) ; 

phiR = -data.anglesBodyFrame(:,PHIR) ;
phiL = data.anglesBodyFrame(:,PHIL) ;
%bodyPitch = data.anglesLabFrame(:,BETA) ;

% -------------------------------------------------------
%% match up fwd/back flip times
[fwdFlipTimesR, fwdFlipTimesL] = ...
    alignFlipTimes(fwdFlipTimesR, fwdFlipTimesL, noMatchCost) ; 
[backFlipTimesR, backFlipTimesL] = ...
    alignFlipTimes(backFlipTimesR, backFlipTimesL, noMatchCost) ; 
%--------------------------------------------------------
%% low-pass butterworth filter for body pitch data
% smoothed body euler angles
[bodyPitch, bodyYaw, bodyRoll] = smoothBodyAngles(data , false, ...
    smoothingParams.pitch_filt_lvl, smoothingParams.yaw_filt_lvl, ...
    smoothingParams.roll_filt_lvl) ;

% combine in matrix and convert to radians
bodyYPR = (pi/180).*[bodyYaw, bodyPitch, bodyRoll] ; 

% get derivatives of euler angles w.r.t. time (vel, accel)
[angleVel, angleAccel] = diffBodyEulerAngles(bodyYPR, dt) ;

c_pitch = fit(t', (180/pi).*bodyYPR(:,2), 'cubicinterp') ; 
c_vel = fit(t', (180/pi).*angleVel(:,2), 'cubicinterp') ; 
c_accel = fit(t', (180/pi).*angleAccel(:,2), 'cubicinterp') ; 
%pitch_filt = filterEulerAngle(bodyPitch,smoothingParams.pitch_filt_lvl) ;
%c_pitch = fit(t',pitch_filt,'cubicinterp');
%pitch_smoothed = c_pitch(t) ;
%pitchVel_smoothed = differentiate(c_pitch, t) ; 
%c_pitchVel = fit(t',pitchVel_smoothed,'linearinterp');

%------------------------------------------
%% spline smooth for stroke angle

sp_phiR = mySplineSmooth(t(~isnan(phiR)),phiR(~isnan(phiR)),phiEstErr) ;
sp_phiL = mySplineSmooth(t(~isnan(phiL)),phiL(~isnan(phiL)),phiEstErr) ;

% -------------------------------------------------------------
%% find front stroke angle corresponding to steady pitch angle
if steadyPhiFrontFlag
    try
        % calculate body pitch acceleration
        %[~, pitch_accel] = differentiate(c_pitch,t) ;
        pitch_accel = c_accel(t) ; 
        
        % get average pitch acceleration for each wingstroke (binned between
        % back flip times)
        meanBackFlipTimes = (backFlipTimesR + backFlipTimesL) ./ 2 ;
        binTimes = (meanBackFlipTimes(2:end) + meanBackFlipTimes(1:end-1))./2 ;
        N_bins = length(meanBackFlipTimes)-1 ;
        
        % also get fwd flip positions
        phiFrontTempR = fnval(sp_phiR,fwdFlipTimesR) ;
        phiFrontTempL = fnval(sp_phiL,fwdFlipTimesL) ;
        
        % loop through bins
        binnedPitchAccel = nan(N_bins,1) ;
        for bin_num = 1:N_bins
            [~, ind1] = min(abs(t - meanBackFlipTimes(bin_num))) ;
            [~, ind2] = min(abs(t - meanBackFlipTimes(bin_num+1))) ;
            binnedPitchAccel(bin_num) = nanmean(pitch_accel(ind1:ind2)) ;
        end
        
        % find wingbeat with least pitch acceleration
        if fitSteadyFlag
            % this version tries to fit relationship between phi front and 
            % pitch  accel (steady value is when linear fit = 0)
            
            % find pert time to exclude from fit
            pert_idx = (binTimes >= 0) & (binTimes <= pulseDuration) ; 
            phiFrontTempMean = (phiFrontTempR + phiFrontTempL)./2 ;
            [~, ~, fwdInd, ~] = alignFlipTimes(fwdFlipTimesR, binTimes) ;
            [c_phi_accel, gof, output] = ...
                fit(phiFrontTempMean(fwdInd & ~pert_idx), ...
                binnedPitchAccel(~pert_idx),'poly1','Robust','LAR') ;

            phiFrontSteady = (-1*c_phi_accel.p2 / c_phi_accel.p1) ;
            
            testDeltaPhiFront = phiFrontTempMean - phiFrontSteady ;
            temp_idx = (fwdFlipTimesR >= manualCorrRange(1)) & ...
                (fwdFlipTimesR <= 0);
            testDeltaPhiFrontRMS = ...
                sqrt(mean((testDeltaPhiFront(temp_idx)).^2)) ;
            if ((testDeltaPhiFrontRMS >= 3) && (gof.dfe < 15)) || ...
                    (testDeltaPhiFrontRMS >= 10) || (gof.rsquare < 0) || ...
                    (output.exitflag < 1)
                phiFrontSteady = nan ;
            end
        else
            % this version finds the wingbeat (prior to pert) with lowest
            % pitch acceleration
            [~, sortIdx] = sort(abs(binnedPitchAccel),'ascend') ;
            binTimesSorted = binTimes(sortIdx) ;
            steadyBinNum = find(binTimesSorted < 0, 1, 'first') ;
            steadyBinTime = binTimesSorted(steadyBinNum) ;
            % find the forward stroke angle that this corresponds to
            
            [~, steadyTimeInd] = min(abs(t - steadyBinTime)) ;
            steady_t = t(steadyTimeInd) ;
            [~, steadyIndR] = min(abs(fwdFlipTimesR - steady_t)) ;
            [~, steadyIndL] = min(abs(fwdFlipTimesL - steady_t)) ;
            phiFrontSteady = (phiFrontTempR(steadyIndR) + ...
                phiFrontTempL(steadyIndL))./2  ;
        end
        
        
        %
        % debug?
        if debugFitFlag
            figure ;
            hold on
            if fitSteadyFlag
                plot(phiFrontTempR(fwdInd), binnedPitchAccel,'ro')
                plot(phiFrontTempL(fwdInd), binnedPitchAccel,'bo')
                plot(phiFrontTempMean(fwdInd), binnedPitchAccel,'ko')
                x_temp = linspace(min(phiFrontTempMean), max(phiFrontTempMean), 20) ;
                plot(x_temp, c_phi_accel(x_temp),'k-')
                xlabel('\Delta \phi_{front} (deg)')
                ylabel('Pitch Accel (deg/s^2)')
                title('Steady state \phi_{front} fit')
            else
                plot(t, pitch_accel)
                axis tight
                plot(binTimes, binnedPitchAccel,'ko-')
                plot(steadyBinTime, binnedPitchAccel(sortIdx(steadyBinNum)),'rx')
                
                xlabel('Time (s)')
                ylabel('Pitch Accel (deg/s^2)')
                title('Wingbeat-averaged pitch accel')
            end
            keyboard
        end
    catch
        phiFrontSteady = nan ;
    end
else
    steady_t = nan ; 
end

%--------------------------------------------
%%  restrict time range
correctedIndR = (fwdFlipTimesR >= manualCorrRange(1)) & ...
    (fwdFlipTimesR <= manualCorrRange(2)) ;
correctedIndL = (fwdFlipTimesL >= manualCorrRange(1)) & ...
    (fwdFlipTimesL <= manualCorrRange(2)) ;
correctedIndBackR = (backFlipTimesR >= manualCorrRange(1)) & ...
    (backFlipTimesR <= manualCorrRange(2)) ;
correctedIndBackL = (backFlipTimesL >= manualCorrRange(1)) & ...
    (backFlipTimesL <= manualCorrRange(2)) ;

correctedInd = correctedIndR | correctedIndL ; 
correctedIndBack = correctedIndBackR | correctedIndBackL ; 

fwdFlipTimesR = fwdFlipTimesR(correctedInd) ;
fwdFlipTimesL = fwdFlipTimesL(correctedInd) ;
backFlipTimesR = backFlipTimesR(correctedIndBack) ;
backFlipTimesL = backFlipTimesL(correctedIndBack) ;

%------------------------------------------
%% if only using one wing...
if ~isfield(data,'oneWing')
    phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
    phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;
elseif strcmp(data.oneWing,'L')
    fwdFlipTimesR = fwdFlipTimesL ;
    phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;
    phiFrontR = phiFrontL ;
elseif strcmp(data.oneWing,'R')
    fwdFlipTimesL = fwdFlipTimesR ;
    phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
    phiFrontL = phiFrontR ;
else
    phiFrontR = fnval(sp_phiR,fwdFlipTimesR) ;
    phiFrontL = fnval(sp_phiL,fwdFlipTimesL) ;
end

%--------------------------------------------------------------------------
%% get mean change in phi front and mean flip times
phiFront = (phiFrontR + phiFrontL) ./ 2 ;
fwdFlipTimes = (fwdFlipTimesR + fwdFlipTimesL ) ./ 2 ;
backFlipTimes = (backFlipTimesR + backFlipTimesL) ./ 2 ; 

if steadyPhiFrontFlag && ~isnan(phiFrontSteady)
    phiFrontAvg = phiFrontSteady ;
    %disp('worked!')
else
    prePertFlipInd = ( fwdFlipTimes < 0 ) ;
    phiFrontAvg = mean(phiFront(prePertFlipInd)) ;
end
deltaPhiFront = (phiFront - phiFrontAvg)' ;
% prePertFlipInd = ( fwdFlipTimes < 0 ) ;
% deltaPhiFront = (((phiFrontR - mean(phiFrontR(prePertFlipInd))) + ...
%     (phiFrontL - mean(phiFrontL(prePertFlipInd))))./2)' ; 
%--------------------------------------------------------------------------
%% make plot to test data?
if debugFlag
    t_min = min([backFlipTimes(1) fwdFlipTimes(1)]) ;
    t_max = max([backFlipTimes(end) fwdFlipTimes(end)]) ;
    t_range = linspace(t_min, t_max, 300) ; 
    
    pitch_smoothed = c_pitch(t_range) ;
    pitchVel = c_vel(t_range) ; 
    pitchAccel = c_accel(t_range) ; 
    %[pitchVel, pitchAccel] = differentiate(c_pitch, t_range) ;
    
    figure ;
    subplot(2,2,1)
    hold on
    plot(1000*t_range, pitch_smoothed,'b-','LineWidth',1.5)
    plot(1000*t, bodyPitch, 'k.')
    %plot(1000*t_peak, pitch_smoothed(ind_peak), 'ro', 'MarkerFaceColor','r')
    set(gca,'xlim', 1000*[t_min t_max])
    ylabel('\theta [deg]')
    
    subplot(2,2,3)
    plot(1000*t_range, pitchVel,'m-','LineWidth',1.5)
    ylabel('Pitch Vel [deg/s]')
    set(gca,'xlim', 1000*[t_min t_max])
    
    subplot(2,2,2)
    plot(1000*t_range, pitchAccel,'r-','LineWidth',1.5)
    hold on
    %plot(1000*t_peak, pitchAccel(ind_peak),'bx')
    set(gca,'xlim',1000*[t_min t_max])
    ylabel('Pitch Accel [deg/s^2]')
    xlabel('Time [ms]')
    
    subplot(2,2,4)
    hold on
    %plot(1000*t_range, c_phiR(t_range), 'r-')
    %plot(1000*t_range, c_phiL(t_range), 'b-')
    plot(1000*fwdFlipTimes, phiFront,'ko-','MarkerFaceColor', 'k')
    %plot(1000*backFlipTimes, phiBack,'k^','MarkerFaceColor', 'k')
    set(gca,'xlim', 1000*[t_min t_max])
    ylabel('\phi [deg]')
end

end