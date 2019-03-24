function [F_lift, F_drag, F_rot, wingTip_smooth, U_t_projHat, spanHat, chordHat_proj, t] = ...
    quasiSteadyForceSam(data,wingSide,plotFlag)
%Sam's attempt at a force calculation going from Sane+Dickinson (2002)
% 
%Translational forces are calculated as:
%
%F_lift = .5 * ( airDensity ) * S * ( U_t )^2 * ( r22_S ) * C_L(alpha)
%F_drag = .5 * ( airDensity ) * S * ( U_t )^2 * ( r22_S ) * C_D(alpha)
%
%   Where:
%       S = projected wing area 
%       U_t = wing tip velocity 
%       r22_S = nondimensional second moment of wing area (Ellington 1984b). Value taken from Sane+Dickingson (2002)
%       C_L, C_D = lift and drag coefficients, respectively. They are functions
%           of angle of attack
%
%Rotational force is calculated as:
%
%F_rot = (C_rot)*(airDensity)*(chord)^2*(span/2)*(r11_v)*(v_hat)*(omega)*U_t
%
%   Where:
%       C_rot = coefficient of rotational force
%       chord = wing mean chord length
%       span = wing span length
%       r11_v = first area moment of virtual mass (Ellington 1984b)
%       v_hat = dimensionless virtual mass (Ellington 1984b)
%       omega = angular velocity of angle of attack (alpha)
%
%Other Outputs:
%   wingTip_smooth = smoothed wing tip positions (voxels)
%   U_t_projHat = direction of wing tip velocity projected in the plane orthogonal to span
%   spanHat = smoothed, normalized span vectors
%   chordHat_proj = smoothed, normalized chord vectors orthogonal to span
%   t = time (s)
%
%Inputs:
%   exprNum = experiment number (7-15 for now)
%   movNum = movie number
%   wingSide = which wing to calculate for. Either 'R' or 'L'
%   plotFlag = 1 if you want plots, 0 if you don't
%
%Example usage:
%   [F_lift, F_drag, F_rot, wingTip_smooth, spanHat, chordHat_proj, t] = quasiSteadyForceSam(7,8,'R',1)
%
%TO DO:
%   -Resolve rotational force into components
%   -Check smoothing (particularly vector component smoothing)
%   -Determine importance of rotational force
%
%THINGS CHANGED WHILE TESTING IDEALIZED WINGSTROKE:
%   -changed WingTipEstErr (line ~156)
%   -removed unit conversion factor from U_t (~line 166)
%
%--------------------------------------------------------------------------

%% Define constants

airDensity = 1.2041 ; %kg/m^3
span = .0025 ; %in meters (was .002)
chord = .0007 ; %in meters
S = pi*(span/2)*(chord/2) ; %area of ellipse
r22_S = .313 ; %Sane+Dickinson (2002) gives .4, Cheng et al. gives .313
r2d = 180/pi ; %radians to degrees conversion
d2r = pi/180 ; %degrees to radians conversion
%wingSide = 'R' ; %'R'

%% Load data
defineConstantsScript
%{
%exprNum = 7 ;
%movNum  = 8 ;

if (movNum<10)
    zstr = '00' ;
elseif (movNum<100)
    zstr = '0' ;
else
    zstr = '';
end

datapath = ['F:\luca\Analysis\pitch down\Expr_' num2str(exprNum) '_mov_' zstr num2str(movNum) '\' ] ;

datafilename = [ datapath ...
    'Expr' num2str(exprNum) 'mov' zstr num2str(movNum) '_Data_manually_corrected.mat' ] ; %

load(datafilename) ;
%}

if (isfield(data,'ignoreFrames'))
    ignoreFrames = data.ignoreFrames ;
else
    ignoreFrames = [] ;
end

%{
if (isfield(data,'correctionTime'))
    correctionTime=data.correctionTime;
elseif(isfield(data,'manualCorrRangeMS'))
    correctionTime=data.manualCorrRangeMS;
end
%}

startTime = data.params.startTrackingTime ;
endTime = data.params.endTrackingTime ;
%pulseLength = data.params.pulseLengthMS/1000 ;
if wingSide == 'R'
    fwdFlipTimes = data.fwdFlipTimesR ; 
    backFlipTimes = data.backFlipTimesR ; 
    %wingPitch = data.anglesLabFrame(:,ETAR) ;
    %wingStroke = data.anglesBodyFrame(:,PHIR) ;
elseif wingSide == 'L'
    fwdFlipTimes = data.fwdFlipTimesL ; 
    backFlipTimes = data.backFlipTimesL ;
    %wingPitch = data.anglesLabFrame(:,ETAL) ;
    %wingStroke = data.anglesBodyFrame(:,PHIL) ;
else 
    disp('Need to fix wingSide')
    return
end
patchColor = [1 1 1 ] * 0.8; %for plotting wingstroke background
fps = data.params.fps ; 
dt = 1/fps ; 


%{
%% Calculate angles (going to use wing pitch as angle of attack)

PsiEstErr = 1 ; 

dt = 1/8000 ;

[anglesLabFrame, anglesBodyFrame, t, newEtaLab, newEtaBody,sp_rho, smoothed_rho, rho_t, rho_samp] ...
    = calcAngles_quick_and_dirty_mk2(data, rollEstErr, false) ;

%Smooth right wing angles in body frame
[anglesMat_Body, smooth_anglesMat_Body, sp_phi_Body, sp_theta_Body, sp_psi_Body ] = ... 
    smoothWingAngles(anglesBodyFrame, t, phiEstErr, thetaEstErr, psiEstErr, wingSide, ignoreFrames);

alpha = fnval(sp_psi_Body,t)' ;
%}
%% Find wing tip velocity

if wingSide == 'R'
    wingTip = data.rightWingTips ;
elseif wingSide == 'L'
    wingTip = data.leftWingTips ;
else
    disp('Need to adjust wingSide')
    return
end

%Need to make sure not to use bad frames. ind are the good indices
t = (startTime:endTime)/fps ;
if (~isempty(ignoreFrames))
    wingTip(ignoreFrames,:) = NaN ; 
end

ind = find(~isnan(wingTip(:,1))) ;
currtvec = t(ind) ;
currwingTip = wingTip(ind,:) ; 

WingTipEstErr = 1.5 ; %1e-6 for idealized wingstroke, 1.5 for real data

[sp_wingTip_X, ~, ~] = mySplineSmooth(currtvec, currwingTip(:,1), WingTipEstErr) ;
[sp_wingTip_Y, ~, ~] = mySplineSmooth(currtvec, currwingTip(:,2), WingTipEstErr) ;
[sp_wingTip_Z, ~, ~] = mySplineSmooth(currtvec, currwingTip(:,3), WingTipEstErr) ;

%evaluate the derivatives of the splines. there's probably a better way to this
wingTip_smooth = [fnval(sp_wingTip_X,t)',fnval(sp_wingTip_Y,t)',fnval(sp_wingTip_Z,t)'] ;
vWingTip = [fnval( fnder(sp_wingTip_X,1), t)', fnval( fnder(sp_wingTip_Y,1), t)', fnval( fnder(sp_wingTip_Z,1), t)'] ;
U_t = vWingTip * 50e-6  ; %convert to meters (now m/s).  

%{
%% effects from rotational velocity
if isfield(data,'thetaDot')
    thetaDot = data.thetaDot ;
else
    thetaDot = 0 ;
end
U_t = U_t + thetaDot*cross(wingTip_smooth*50e-6,repmat([1 0 0],size(wingTip_smooth,1),1));
%}
%{
%% Body velocity
%This is done to make sure forces are calculated in body frame (Dickinson's
%stage doesn't move
CMEstErr = 1 ;
CM = data.bodyCM ;

tol = CMEstErr ;
[sp_CM_X, CM_X_smooth, ~] = mySplineSmooth(t, CM(:,1), tol) ;
[sp_CM_Y, CM_Y_smooth, ~] = mySplineSmooth(t, CM(:,2), tol) ;
[sp_CM_Z, CM_Z_smooth, ~] = mySplineSmooth(t, CM(:,3), tol) ;

%evaluate the derivatives of the splines. there's probably a better way to this
vBody = [fnval( fnder(sp_CM_X,1), t)', fnval( fnder(sp_CM_Y,1), t)', fnval( fnder(sp_CM_Z,1), t)'] ;
vBody = vBody * 50e-6 ; %convert to meters (now m/s)

U_t_Body = U_t - vBody ; 
%}
%% Find alpha by taking angle between velocity and chord

if wingSide == 'R'
    chordHat = data.rightChordHats ;
    spanHat = data.rightSpanHats ;
elseif wingSide == 'L'
    chordHat = data.leftChordHats ;
    spanHat = data.leftSpanHats ;
else
    disp('Need to adjust wingSide')
    return
end

%smooth span vectors
currSpanHat = spanHat(ind,:) ;
currSpanTheta = 90*ones(size(currSpanHat(:,3))) - (r2d)*asin(currSpanHat(:,3)) ; %these are standard spherical polar coordinates now
currSpanPhi = (r2d)*unwrap(atan2(currSpanHat(:,2),currSpanHat(:,1))) ;
spanThetaEstErr = .5 ; %degrees
spanPhiEstErr = 1 ; %degrees

[sp_spanTheta, ~, ~] = mySplineSmooth(currtvec, currSpanTheta, spanThetaEstErr) ;
[sp_spanPhi, ~, ~] = mySplineSmooth(currtvec, currSpanPhi, spanPhiEstErr) ;

spanHatTheta = fnval(sp_spanTheta,t)' ;
spanHatPhi = fnval(sp_spanPhi,t)' ;
spanHat = [cos(d2r*spanHatPhi).*sin(d2r*spanHatTheta), ...
    sin(d2r*spanHatPhi).*sin(d2r*spanHatTheta), cos(d2r*spanHatTheta)] ;

%Remove spanwise component of velocity
U_t_proj = U_t ;%- repmat(dot(U_t,spanHat,2),1,3).*spanHat ;
U_t_projNorm = sqrt(sum(U_t_proj.*U_t_proj,2)) ;
U_t_projHat = U_t_proj ./ repmat(U_t_projNorm,1,3) ;

%smooth chord vectors
chordHat_proj = chordHat ; %- repmat(dot(chordHat,spanHat,2),1,3).*spanHat ;
currChordHat_proj = chordHat_proj(ind,:) ; 
currChordHat_projTheta = 90*ones(size(currChordHat_proj(:,3))) - (r2d)*unwrap(asin(currChordHat_proj(:,3))) ;
currChordHat_projPhi = (r2d)*unwrap(atan2(currChordHat_proj(:,2),currChordHat_proj(:,1))) ;
chordThetaEstErr = 1 ;
chordPhiEstErr = 2 ;

[sp_chordTheta, ~, ~] = mySplineSmooth(currtvec, currChordHat_projTheta, chordThetaEstErr) ;
[sp_chordPhi, ~, ~] = mySplineSmooth(currtvec, currChordHat_projPhi, chordPhiEstErr) ;

chordHat_projTheta = fnval(sp_chordTheta,t)' ; 
chordHat_projPhi = fnval(sp_chordPhi,t)' ; 
chordHat_proj = [cos(d2r*chordHat_projPhi).*sin(d2r*chordHat_projTheta), ...
    sin(d2r*chordHat_projPhi).*sin(d2r*chordHat_projTheta), cos(d2r*chordHat_projTheta)] ;


alpha = (180/pi)*acos(sum(U_t_projHat .* chordHat_proj,2)) ; 

%Plot alpha
if plotFlag
    halpha = figure ;
    
    set(gca,'fontsize',14);
    plot(t*1000,alpha,'LineWidth',2,'Color',[.7 0 0]) ;
    hold on ;
    box on; grid on;
    plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
    axis tight
    xlabel('Time [ms]')
    ylabel('\alpha [deg]')
    title('Angle of Attack')
end
%% Calculate force magnitudes

N = length(alpha) ; 
C_lift = 0.225*ones(N,1) + 1.58*sin((pi/180)*(2.13*alpha - 7.2*ones(N,1))) ;
C_drag = 1.92*ones(N,1) - 1.55*cos((pi/180)*(2.04*alpha - 9.82*ones(N,1))) ; 

F_lift = .5 * ( airDensity ) * S * ( r22_S ) * C_lift .* U_t_projNorm.^2 ;
F_drag = .5 * ( airDensity ) * S  * ( r22_S )* C_drag .* U_t_projNorm.^2 ;

%% Make figures
if plotFlag
    hforce = figure('position',[140 50 1100 900]) ;
    
    subplot(3,1,1);
    %hold on;
    set(gca,'fontsize',14) ;
    plot(t*1000, F_lift, 'b', 'LineWidth', 2)
    grid on ; box on;
    hold on ;
    plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
    ylabel('Lift Force [N]')
    title(['Translation Force Calculation ' wingSide])
    axis tight
    
    subplot(3,1,2)
    set(gca,'fontsize',14) ;
    plot(t*1000, F_drag, 'r', 'LineWidth', 2)
    grid on ; box on;
    hold on ;
    plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
    ylabel('Drag Force [N]')
    axis tight
    
    subplot(3,1,3)
    set(gca,'fontsize',14) ;
    plot(t*1000, sqrt(F_drag.^2+F_lift.^2), 'Color',[.5 0 .5], 'LineWidth', 2)
    grid on ; box on;
    hold on ;
    plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
    ylabel('Total Force [N]')
    axis tight
    xlabel('Time [ms]')
    
    hvelocity = figure('position',[140 50 1100 900]) ;
    
    subplot(3,1,1)
    set(gca,'fontsize',14);
    plot(t*1000,U_t(:,1),'r','linewidth',2);
    grid on ; box on;
    hold on ;
    plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
    title(['Wing Tip Velocities ' wingSide])
    ylabel('v_x [m/s]')
    axis tight
    
    subplot(3,1,2)
    set(gca,'fontsize',14);
    plot(t*1000,U_t(:,2),'b','linewidth',2);
    grid on ; box on;
    hold on ;
    plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
    ylabel('v_y [m/s]')
    axis tight
    
    subplot(3,1,3)
    set(gca,'fontsize',14);
    plot(t*1000,U_t(:,3),'Color',[0 .5 0],'linewidth',2);
    grid on ; box on;
    hold on ;
    plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
    ylabel('v_z [m/s]')
    xlabel('Time [ms]')
    axis tight
end 
%% Calculate the rotational force

%Define some of the constants
x_0 = .25 ;                 %rotational axis for wing (see Sane+Dickinson 2002)
C_rot = pi*(.75 - x_0) ;
v_hat = 1.1;                %see Ellington 1984b
r11_v = 1.72*(r22_S)^1.07;  %see Ellington 1984b

%Get the angular velocity ( d(alpha)/dt )
%wingPitchEstErr = 4 ;
%[sp_wingPitch, ~, ~] = mySplineSmooth(t, wingPitch, wingPitchEstErr) ;
%omega = fnval( fnder(sp_wingPitch,1), t)' ; 

%wingStrokeEstErr = 1 ;
%[sp_wingStroke, ~, ~] = mySplineSmooth(t, wingStroke, wingStrokeEstErr) ;
%phidot = -fnval( fnder(sp_wingStroke,1), t)' ; 


%Throw it all together
F_rot = zeros(size(F_lift)) ;%(C_rot)*(airDensity)*(chord)^2*(span/2)*(r11_v)*(v_hat)*...
    %(pi/180)*omega.*U_t_projNorm.*sign(phidot) ; %Eq. 11 in Sane+Dickinson 11 w/ modifications from Ellington

%Plot rotational force
if plotFlag
    hrotforce = figure('position',[140 50 1100 300]) ; 
    
    set(gca,'fontsize',14);
    plot(t*1000,F_rot,'Color', [0 .5 0], 'linewidth',2);
    hold on ;
    grid on ; box on;
    plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
    title(['Rotational Force ' wingSide])
    ylabel('Force [N]')
    xlabel('Time [ms]')
    axis tight
end

%% Misc. functions
%{
function C_L = lift_coeff(alpha)
    C_L = 0.225 + 1.58*sind(2.13*alpha - 7.2) ;
end
function C_D = drag_coeff(alpha)
    C_D = 1.92 - 1.55*cosd(2.04*alpha - 9.82) ;
end
%}