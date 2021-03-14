function [anglesMat, smooth_anglesMat, sp_phi, sp_theta, sp_psi ] = ...
    smoothWingAngles(data, wingSide, frameType)
%--------------------------------------------------------------------------
%Takes the output of calculate angle and returns arrays of angles, smoothed
%angles, and splines. Order is [phi, theta, psi]
%
%wingSide is either 'L' or 'R'
%
%--------------------------------------------------------------------------
%% params
if ~exist('frameType','var') || isempty(frameType)
    frameType = 'Body' ;
end
debugFlag = false ;
debugFlag2 = false ;

defineConstantsScript
smoothingParams = setSmoothingParams() ;

% spline fit
phiEstErr = smoothingParams.phi_est_err ;
thetaEstErr = smoothingParams.theta_est_err ;
psiEstErr = smoothingParams.eta_est_err ;
%psiEstErr2 = 1 ;

% low pass filter (not currently used)
phi_filt_lvl = smoothingParams.phi_filt_lvl ;
theta_filt_lvl = smoothingParams.theta_filt_lvl ;
eta_filt_lvl = smoothingParams.eta_filt_lvl ;

% savitzky golay
k_phi = smoothingParams.phi_sgolay_k ;
k_theta = smoothingParams.theta_sgolay_k ;
k_psi = smoothingParams.eta_sgolay_k ;
f_phi = smoothingParams.phi_sgolay_f ;
f_theta = smoothingParams.theta_sgolay_f ;
f_psi = smoothingParams.eta_sgolay_f ;

% hampel filter
hampel_k = smoothingParams.wing_hampel_k ;
hampel_nsigma = smoothingParams.wing_hampel_nsigma ;
%--------------------------------------------------------------------------
%% load data
dt = 1/8000 ;
angles = data.(['angles' frameType 'Frame']) ;
t = (data.params.startTrackingTime : data.params.endTrackingTime )*dt ;

if isfield(data,'manualCorrRangeMS')
    tlim = data.manualCorrRangeMS / 1000 ;
else
    tlim = [t(1) t(end)] ;
end

if isfield(data,'ignoreFrames')
    ignoreFrames = data.ignoreFrames ;
else
    ignoreFrames = [] ;
end

switch wingSide
    case 'L'
        phi = angles(:, PHIL) ;
        theta = angles(:, THETAL) ;
        psi = (180/pi)*unwrap((pi/180)*angles(:, ETAL)) ;
        ignoreIndPhi = unique([find(isnan(angles(:, PHIL))==1)' ignoreFrames]) ;
    case 'R'
        phi = angles(:,PHIR) ;
        if ~strcmp(frameType,'Lab')
            phi = -1.*phi ;
        end
        theta = angles(:, THETAR) ;
        psi = (180/pi)*unwrap((pi/180)*angles(:, ETAR)) ;
        ignoreIndPhi = unique([find(isnan(angles(:, PHIR))==1)' ignoreFrames]) ;
    otherwise
        disp('Need to specify wing side')
        return
end

%check wing pitch to see if unwrap worked
for pp = 1:length(psi)
    while psi(pp) < 0
        psi(pp) = psi(pp) + 360 ;
    end
    while psi(pp) > 360
        psi(pp) = psi(pp) - 360 ;
    end
end

%try to detect/remove outliers
[~, hampelTheta] = hampel(theta, hampel_k, hampel_nsigma) ;
hampelThetaInd = find(hampelTheta) ;
[~, hampelPsi] = hampel(psi, hampel_k, hampel_nsigma) ;
hampelPsiInd = find(hampelPsi) ;

ignoreIndPhi = unique([ignoreIndPhi, hampelThetaInd', hampelPsiInd' ]) ;
ignoreIndPsi = unique([ignoreFrames, hampelPsiInd' ]) ;
if (~isempty(ignoreIndPhi))
    phi(ignoreIndPhi) = NaN ;
    theta(ignoreFrames) = NaN ;
    psi(ignoreIndPsi) = NaN ;
end


%should already be done by quick_and_dirty
%[~, ~, ~, ~, ~, ~, badIndicesR] = findWingFlipTimes_mk3 (t, phi, false);
%phi(badIndicesR) = NaN ;
%theta(badIndicesR) = NaN ;
%psi(badIndicesR) = NaN ;

%% interpolate ignore frames
nanIndPhi = isnan(phi) ;
notNanIndPhi = ~isnan(phi) ;
c_phi = fit(t(notNanIndPhi)',  phi(notNanIndPhi),'smoothingspline') ;
% c_phi = fit(t(notNanIndPhi)',  phi(notNanIndPhi),'cubicinterp') ;
phi_interp = c_phi(t(nanIndPhi)) ;
%phi_interp = interp1(t(notNanIndPhi), phi(notNanIndPhi), t(nanIndPhi),'spline') ;

nanIndTheta = isnan(theta) ;
notNanIndTheta = ~isnan(theta) ;
c_theta = fit(t(notNanIndTheta)',  theta(notNanIndTheta),'smoothingspline') ;
% c_theta = fit(t(notNanIndTheta)',  theta(notNanIndTheta),'cubicinterp') ;
theta_interp = c_theta(t(nanIndTheta)) ;
%theta_interp = interp1(t(notNanIndTheta), theta(notNanIndTheta), t(nanIndTheta),'spline') ;

nanIndPsi = isnan(psi) ;
notNanIndPsi = ~isnan(psi) ;
c_psi = fit(t(notNanIndPsi)',  psi(notNanIndPsi),'smoothingspline') ;
% c_psi = fit(t(notNanIndPsi)',  psi(notNanIndPsi),'cubicinterp') ;
psi_interp = c_psi(t(nanIndPsi)) ;
%psi_interp = interp1(t(notNanIndPsi), psi(notNanIndPsi), t(nanIndPsi),'spline') ;

if debugFlag
    
    t1Ind = find(t == tlim(1)) ;
    t2Ind = find(t == tlim(2)) ;
    
    figure ;
    subplot(3,1,1)
    hold on
    plot(1000*t(notNanIndPhi), phi(notNanIndPhi), 'k.')
    plot(1000*t(nanIndPhi), phi_interp, 'rx')
    set(gca,'xlim',tlim*1000)
    set(gca,'ylim', [min(phi(t1Ind:t2Ind))-2 max(phi(t1Ind:t2Ind)+2)])
    title('\phi raw data + interpolated values')
    ylabel('\phi [deg]')
    
    subplot(3,1,2)
    hold on
    plot(1000*t(notNanIndTheta), theta(notNanIndTheta), 'k.')
    plot(1000*t(nanIndTheta), theta_interp, 'rx')
    set(gca,'xlim',tlim*1000)
    set(gca,'ylim', [min(theta(t1Ind:t2Ind))-2 max(theta(t1Ind:t2Ind)+2)])
    %axis tight
    title('\theta raw data + interpolated values')
    ylabel('\theta [deg]')
    
    subplot(3,1,3)
    hold on
    plot(1000*t(notNanIndPsi), psi(notNanIndPsi), 'k.')
    plot(1000*t(nanIndPsi), psi_interp, 'rx')
    set(gca,'xlim',tlim*1000)
    set(gca,'ylim', [min(psi(t1Ind:t2Ind))-2 max(psi(t1Ind:t2Ind)+2)])
    %axis tight
    title('\psi raw data + interpolated values')
    xlabel('Time [ms]')
    ylabel('\psi [deg]')
end

phi_comb = phi ; phi_comb(nanIndPhi) = phi_interp ;
theta_comb = theta ; theta_comb(nanIndTheta) = theta_interp ;
psi_comb = psi ; psi_comb(nanIndPsi) = psi_interp ;

anglesMat = [phi_comb,theta_comb,psi_comb] ;

%% smooth angles
% phi_smooth = sgolayfilt(phi_comb, k_phi, f_phi) ;
% theta_smooth = sgolayfilt(theta_comb, k_theta, f_theta) ;
% psi_smooth = sgolayfilt(psi_comb, k_psi, f_psi) ;

phi_smooth = smooth(phi_comb, f_phi, 'sgolay', k_phi) ;
theta_smooth = smooth(theta_comb, f_theta, 'sgolay', k_theta) ;
psi_smooth = smooth(psi_comb, f_psi, 'sgolay', k_psi) ;



phi_filt = filterEulerAngle(phi_comb,phi_filt_lvl) ;
theta_filt = filterEulerAngle(theta_comb, theta_filt_lvl) ;
psi_filt = filterEulerAngle(psi_comb, eta_filt_lvl) ;

[sp_phi, ~, ~] =  mySplineSmooth(t, phi_comb, phiEstErr) ;
[sp_theta, ~, ~] =  mySplineSmooth(t, theta_comb, thetaEstErr) ;
[sp_psi, ~, ~] =  mySplineSmooth(t, psi_smooth, psiEstErr) ;

if debugFlag2
    
    t1Ind = find(t == tlim(1)) ;
    t2Ind = find(t == tlim(2)) ;
    
    %under construction
    phi_smooth2 = smooth(phi, 9, 'rloess') ;
    theta_smooth2 = smooth(theta_comb, 7, 'rloess') ;
    psi_smooth2 = smooth(psi_comb, 7, 'rloess') ;
    
    figure ;
    subplot(3,1,1)
    hold on
    plot(t*1000,phi_comb,'k.')
    plot(t*1000,phi_smooth2,'r-')
    plot(t*1000, phi_smooth, 'b-')
    plot(t*1000, phi_filt,'-','Color',[0 .7 0])
    set(gca,'xlim',tlim*1000)
    set(gca,'ylim', [min(phi(t1Ind:t2Ind))-2 max(phi(t1Ind:t2Ind)+2)])
    title('sgolay tests')
    ylabel('\phi [deg]')
    
    
    subplot(3,1,2)
    hold on
    plot(t*1000,theta_comb,'k.')
    plot(t*1000,theta_smooth2,'r-')
    plot(t*1000, theta_smooth, 'b-')
    %plot(t*1000, medfilt1(theta_comb),'-','Color',[0 .7 0])
    set(gca,'xlim',tlim*1000)
    set(gca,'ylim', [min(theta(t1Ind:t2Ind))-2 max(theta(t1Ind:t2Ind)+2)])
    ylabel('\theta [deg]')
    
    subplot(3,1,3)
    hold on
    plot(t*1000,psi_comb,'k.')
    plot(t*1000,psi_smooth2,'r-')
    plot(t*1000, psi_smooth, 'b-')
    %plot(t*1000, psi_filt,'-','Color',[0 .7 0])
    set(gca,'xlim',tlim*1000)
    set(gca,'ylim', [min(psi(t1Ind:t2Ind))-2 max(psi(t1Ind:t2Ind)+2)])
    %axis tight
    xlabel('Time [ms]')
    ylabel('\psi [deg]')
end

smooth_anglesMat = [phi_smooth'; theta_smooth'; psi_smooth'] ;


end






