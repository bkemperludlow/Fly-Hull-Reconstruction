%--------------------------------------------------------------------------
% function to fetch relevant data for PI controller fitting of yaw
% perturbations
%--------------------------------------------------------------------------
function [deltaAlphaMean, c_yaw, midWingBeatTimes, deltaAlphaSEM] = ...
    getYawControllerData(data, debugFlag)

defineConstantsScript
smoothingParams = setSmoothingParams ;

if isfield(data,'manualCorrRangeMS')
    manualCorrRangeMS_start = max([data.manualCorrRangeMS(1) -10]) ;
    manualCorrRangeMS_end = min([data.manualCorrRangeMS(2) 40]) ;
    manualCorrRangeMS = [manualCorrRangeMS_start, manualCorrRangeMS_end ] ;
else
    manualCorrRangeMS = [-10 60] ;
end
%manualCorrRangeMS = [-10 60] ;
manualCorrRange = manualCorrRangeMS / 1000 ; 

%--------------------------------------------------------------------------
% wing flip times
fwdFlipTimesR = data.fwdFlipTimesR ;
fwdFlipTimesL = data.fwdFlipTimesL ;
backFlipTimesR = data.backFlipTimesR ;
backFlipTimesL = data.backFlipTimesL ;

% wing flip time indices
fwdFlipIndR = data.fwdFlipIndR ;
fwdFlipIndL = data.fwdFlipIndL ;
backFlipIndR = data.backFlipIndR ;
backFlipIndL = data.backFlipIndL ;

% find wing flip times within time range
correctedIndR = find(fwdFlipTimesR > manualCorrRange(1) & fwdFlipTimesR < manualCorrRange(2)) ;
correctedIndL = find(fwdFlipTimesL > manualCorrRange(1) & fwdFlipTimesL < manualCorrRange(2)) ;
correctedIndBackR = find(backFlipTimesR > manualCorrRange(1) & backFlipTimesR < manualCorrRange(2)) ;
correctedIndBackL = find(backFlipTimesL > manualCorrRange(1) & backFlipTimesL < manualCorrRange(2)) ;

fwdFlipTimesR = fwdFlipTimesR(correctedIndR) ;
fwdFlipTimesL = fwdFlipTimesL(correctedIndL) ;
backFlipTimesR = backFlipTimesR(correctedIndBackR) ;
backFlipTimesL = backFlipTimesL(correctedIndBackL) ;

fwdFlipIndR = fwdFlipIndR(correctedIndR) ;
fwdFlipIndL = fwdFlipIndL(correctedIndL) ;
backFlipIndR = backFlipIndR(correctedIndBackR) ;
backFlipIndL = backFlipIndL(correctedIndBackL) ;

%--------------------------------------------------------------------------
% get time and angle variables

t = (data.params.startTrackingTime : data.params.endTrackingTime ) / 8000 ;
bodyYaw = data.anglesLabFrame(:,PSI) ;

% smooth wing angles
[~, smooth_anglesMat_L, ~, ~, sp_etaL ] = ... 
    smoothWingAngles(data, 'L') ; 

[~, smooth_anglesMat_R, ~, ~, sp_etaR ] = ... 
    smoothWingAngles(data, 'R') ; 


etaR = smooth_anglesMat_R(3,:) ;
etaL = smooth_anglesMat_L(3,:) ;

%low-pass butterworth filter for data
yaw_filt = filterEulerAngle(bodyYaw,smoothingParams.yaw_filt_lvl) ;
c_yaw = fit(t',yaw_filt,'cubicinterp');

%--------------------------------------------------------------------------
% get angle of attack difference for right and left wing
% (right - left) for forward stroke, (left-right) for back stroke 

D = pdist2(fwdFlipTimesL', fwdFlipTimesR') ;
[match_ind, ~] = assignmentallpossible(D) ;

numFwdFlipTimes = min([length(fwdFlipTimesR), length(fwdFlipTimesL)]) ; 

deltaAlphaMean = nan(numFwdFlipTimes-1,1) ; 
deltaAlphaSEM = nan(numFwdFlipTimes-1,1) ; 
midWingBeatTimes = nan(numFwdFlipTimes-1,1) ;

DELTA = 2 ; 

for i = 1:(numFwdFlipTimes - 1)
    
    
    backFlipIndL_curr = backFlipIndL(find(backFlipIndL > fwdFlipIndL(match_ind == i),1,'first')) ;
    backFlipIndR_curr = backFlipIndR(find(backFlipIndR > fwdFlipIndR(i),1,'first')) ;
    backFlipInd_curr = round((backFlipIndL_curr + backFlipIndR_curr)/2) ; 
    
    fwdFlipIndL_curr = fwdFlipIndL(match_ind == i) ; 
    fwdFlipIndR_curr = fwdFlipIndR(i) ;  
    fwdFlipInd_curr = round((fwdFlipIndL_curr + fwdFlipIndR_curr)/2) ;
    
    fwdFlipIndL_next = fwdFlipIndL(match_ind==(i+1)) ; 
    fwdFlipIndR_next = fwdFlipIndR(i+1) ; 
    fwdFlipInd_next = round((fwdFlipIndL_next + fwdFlipIndR_next)/2) ;
    

    ind_back = (fwdFlipInd_curr+DELTA):(backFlipInd_curr-DELTA) ; 
    etaL_back = fnval(sp_etaL,t(ind_back)) ;
    etaR_back = fnval(sp_etaR,t(ind_back)) ;
    
    ind_fwd = (backFlipInd_curr+DELTA):(fwdFlipInd_next-DELTA) ;
    etaL_fwd = fnval(sp_etaL,t(ind_fwd)) ;
    etaR_fwd = fnval(sp_etaR,t(ind_fwd)) ;
   
    
    etaL_comb = [etaL_back, -1*etaL_fwd] ;
    etaR_comb = [etaR_back, -1*etaR_fwd] ;
    
    deltaAlpha = etaL_comb - etaR_comb ; 
    deltaAlphaMean(i) = median(deltaAlpha) ; 
    deltaAlphaStd = std(deltaAlpha) ; 
    deltaAlphaSEM(i) = deltaAlphaStd/length(deltaAlpha) ; 
    
    midWingBeatTimes(i) = (t(fwdFlipInd_curr) + t(fwdFlipInd_next))/2 ; 

%     deltaAlphaMean = [deltaAlphaMean ; mean(etaL_back - etaR_back) ; ...
%         mean(etaR_fwd - etaL_fwd)] ; 
%     deltaAlphaSEM = [deltaAlphaSEM ; std(etaL_back - etaR_back)/length(etaL_back) ; ...
%         std(etaR_fwd - etaL_fwd)/length(etaR_fwd)] ; 
%     midWingBeatTimes = [midWingBeatTimes ; ...
%         (t(fwdFlipInd_curr) + t(backFlipInd_curr))/2 ; ...
%         (t(backFlipInd_curr) + t(fwdFlipInd_next))/2 ] ; 
    
end
   

if debugFlag
    t_min = min([backFlipTimesR(1), fwdFlipTimesR(1), ...
        backFlipTimesL(1), fwdFlipTimesL(1)]) ;
    t_max = max([backFlipTimesR(end) fwdFlipTimesR(end),...
        backFlipTimesL(end) fwdFlipTimesL(end)]) ;
    t_range = linspace(t_min, t_max, 300) ; 
    
    yaw_smoothed = c_yaw(t_range) ;
    [yawVel, yawAccel] = differentiate(c_yaw, t_range) ;
    
    figure ;
    subplot(2,2,1)
    hold on
    plot(1000*t_range, yaw_smoothed,'b-','LineWidth',1.5)
    plot(1000*t, bodyYaw, 'k.')
    %plot(1000*t_peak, pitch_smoothed(ind_peak), 'ro', 'MarkerFaceColor','r')
    ylabel('\psi [deg]')
    
    subplot(2,2,3)
    plot(1000*t_range, yawVel,'m-','LineWidth',1.5)
    ylabel('Yaw Vel [deg/s]')
    
    subplot(2,2,2)
    plot(1000*t_range, yawAccel,'r-','LineWidth',1.5)
    hold on
    %plot(1000*t_peak, pitchAccel(ind_peak),'bx')
    %set(gca,'xlim',1000*[t_min t_max])
    ylabel('Yaw Accel [deg/s^2]')
    xlabel('Time [ms]')
    
    subplot(2,2,4)
    hold on
    %plot(1000*t_range, fnval(sp_etaR,t_range), 'r-')
    %plot(1000*t_range, fnval(sp_etaL,t_range), 'b-')
    errorbar(1000*midWingBeatTimes, deltaAlphaMean,deltaAlphaSEM,'ksq',...
        'MarkerFaceColor', 'k')
    %plot(1000*backFlipTimes, phiBack,'k^','MarkerFaceColor', 'k')
    set(gca,'xlim', 1000*[t_min t_max])
    ylabel('\Delta\alpha [deg]')
    
    figure ; 
    hold on
    plot(1000*t_range, fnval(sp_etaR,t_range), 'r-')
    plot(1000*t_range, fnval(sp_etaL,t_range), 'b-')
    plotWingstrokeBackground(gca, 1000*backFlipTimesR, 1000*fwdFlipTimesR, 0.7*[1,1,1], true);
    set(gca,'xlim', 1000*[t_min t_max])
    ylabel('\Delta\alpha [deg]')
end




end