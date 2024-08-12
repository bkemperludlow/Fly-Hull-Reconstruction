%--------------------------------------------------------------------------
% function to calculate body and wing angles from the data structure output
% by hullAnalysis. Based on calcAngles_quick_and_dirty_mk*.m
%
% this version uses a vector-based calculation of the wing pitch angles eta.
% this calculation is done in the function calcEta.m :
% calcEta(s, c, leftRight)
%
%
% each matrix contains 8 columns for
% beta (body), psi (body), phiR, phiL, thetaR, thetaL, etaR, etaL
%
% subfunction used:
% calcBodyEulerAngles
% calcEta
%--------------------------------------------------------------------------
function [anglesLabFrame, anglesBodyFrame, t, newEtaLab, newEtaBody, sp_rho,...
    smoothed_rho, rho_t, rho_samp, rotM_YP, rotM_roll, largePertFlag] = ...
    calcAnglesRaw_Sam(data, plotFlag,largePertFlag)
%--------------------------------------------------------------------------
%% params and inputs
%
if (~exist('plotFlag','var'))
    plotFlag = false ;
end
if (~exist('largePertFlag','var'))
    largePertFlag = false ;
end
% when putting the body in consistent reference frame, make the body axis
% pitched up by 45 degrees from the horizontal plane. this should make the
% stroke plane roughly horizontal in the transformed frame.
%
% NB: should I try to get a better estimate for this? in principle, it
% would affect stroke vs deviation wing angles
thetaB0 = 45 ;
thetaB0rad = thetaB0 * pi / 180 ;

% misc useful params
DEG2RAD = pi / 180 ;
RAD2DEG   = 180 / pi ;
tol = 1e-12 ;
fps = data.params.fps ;

% the angles will ultimately be stored in Nx8 or Nx9 matrices--the script
% below defines variables whose names are the greek letters corresponding
% to the angles, so it's easier to access the right angle
defineConstantsScript

%--------------------------------------------------------------------------
%% initialize data containers
anglesLabFrame  = zeros(data.Nimages, 9);
anglesBodyFrame = zeros(data.Nimages, 8);

newEtaLab  = zeros(data.Nimages,2) ;
newEtaBody = zeros(data.Nimages,2) ;

% timing info
if (isfield(data,'startAnalysisTimeMS'))
    startTime = data.startAnalysisTimeMS * fps / 1000 ;
    endTime = data.endAnalysisTimeMS * fps / 1000 ;
    Np = endTime - startTime + 1 ;
else
    startTime = data.params.startTrackingTime ;
    endTime   = data.params.endTrackingTime ;
    Np = data.Nimages ;
end

t = (startTime:endTime) / fps  ; % in SEC
startFrame = startTime - data.params.startTrackingTime + 1 ;
endFrame   = startFrame + Np - 1 ;
frameIndices = startFrame:endFrame ;


allT   = (0:data.Nimages-1) + data.params.startTrackingTime ;
allT   = allT / data.params.fps ; % in SEC
% -------------------------------------------------------------------------
%% calculate raw psi (yaw) and raw beta (pitch) for all frames (IN DEGREES)
% ----------------------------------------------------------
if largePertFlag
    % in cases where gimbal lock may be a problem, use a frame-by_frame
    % estimate for the body pitch and yaw angles
    [rawBeta, rawPsi, rotM_YP, largePertFlag] = ...
        calcPitchLargePert(data, plotFlag) ;
else
    % otherwise proceed as normal
    rawPsi  = zeros(data.Nimages,1) ; % yaw
    rawBeta = zeros(data.Nimages,1) ; % pitch
    rotM_YP = zeros(3,3,data.Nimages) ; % rotation matrices that will unyaw + unpitch
    for k=1:data.Nimages
        AHat = data.AHat(k,:) ;
        rawPsi(k)  = atan2(AHat(2),AHat(1)); % body angle with respect to x axis 
        rawBeta(k) = asin(AHat(3));  % body angle with respect to the horizon 
        rotM_YP(:,:,k) = eulerRotationMatrix(rawPsi(k), rawBeta(k),0) ; 
    end
    % convert rawPsi and rawBeta to degrees
    rawPsi = RAD2DEG*rawPsi ; 
    rawBeta = RAD2DEG*rawBeta ; 
end
% -------------------------------------------------------------------------
%% calculate body roll
% calc roll angles at data.rhoTimes
if (~isfield(data,'rhoTimes')) || (isempty(data.rhoTimes))
    %[rhoTimes, rollVectors] = estimateRollVectors(data);
    disp('data.rhoTimes is not defined.') ;
    rhoFlag = false ;
else
    rhoTimes    = data.rhoTimes ;
    rollVectors = data.rollVectors ;
    rhoFlag = true ;
end

if rhoFlag && isfield(data, 'newRhoSamp')
    % if possible, uses 'newRhoSamp', estimated in the script for large 
    % perturbation corrections
    [smoothed_rho, sp_rho, rho_t, rho_samp] = calcRollNewRhoSamp(data,t) ; 
elseif rhoFlag && ~isfield(data, 'newRhoSamp')
    % if we have estimates of roll vector, determine body roll angle
    [smoothed_rho, sp_rho, rho_t, rho_samp, rotM_roll] = ...
        calcBodyRoll(rhoTimes, rollVectors, t, rotM_YP, data.params,...
         largePertFlag) ;
else
    smoothed_rho = 0 ; %rho0 ;
    sp_rho = [] ;
    rho_t = allT ;
    rho_samp = 0 ;
end

%--------------------------------------------------------------------------
%% now calculate wing angles.
% to do so, we move everything to the body frame by undoing the yaw, pitch,
% and roll rotations. then calculate wing angles in the new frame. We also
% calculate wing angles in the lab frame
for k=1:data.Nimages
    
    % ------------------
    %% COPY RELEVANT DATA
    % ------------------
    cb = data.bodyCM(k,:)' ;
    cr = data.rightWingCM(k,:)' ;
    cl = data.leftWingCM(k,:)' ;
    
    rightChordHat = data.rightChordHats(k,:)' ;
    leftChordHat  = data.leftChordHats(k,:)' ;
    rightSpanHat  = data.rightSpanHats(k,:)' ;
    leftSpanHat   = data.leftSpanHats(k,:)' ;
    
    % get wing pitch right away
    newEtaLab(k,1) = calcEta(rightSpanHat, rightChordHat,'right') ;
    newEtaLab(k,2) = calcEta(leftSpanHat, leftChordHat,'left') ;
    
    etaRdeg = newEtaLab(k,1)  ;
    etaLdeg = newEtaLab(k,2)  ;
    
    % Use the smoothed version of the roll angle rho
    bodyRollAngle = smoothed_rho(k) ;
    
    % ---------------------------------------------------------------------
    %% CALCULATE ANGLES IN THE LAB FRAME
    % ---------------------------------------------------------------------
    % body
    psiDeg  = rawPsi(k) ;
    betaDeg = rawBeta(k) ;
    
    % right wing
    phiRdeg   = atan2(rightSpanHat(2),rightSpanHat(1)) * RAD2DEG ;
    thetaRdeg = asin(rightSpanHat(3)) * RAD2DEG;
    
    % left wing
    phiLdeg   = atan2(leftSpanHat(2), leftSpanHat(1)) * RAD2DEG;
    thetaLdeg = asin(leftSpanHat(3)) * RAD2DEG;
    
    % store result
    anglesLabFrame(k,:) = ...
        [psiDeg betaDeg phiRdeg thetaRdeg etaRdeg phiLdeg...
        thetaLdeg etaLdeg bodyRollAngle] ;
    
    % -----------------------------------------------
    %% CALCULATE ANGLES IN THE BODY FRAME OF REFERENCE
    % -----------------------------------------------
    
    % first rotation matrix bring (xlab, ylab, zlab) to
    % (xbody ybody zbody) such that xbody is AHat
    
%     % convert to radians
%     phiB   = psiDeg*deg2rad ;
%     thetaB = betaDeg*deg2rad ;
%     psiB   = bodyRollAngle*deg2rad ;
%     
%     % rotation matrix to strict body axis
%     M1 = eulerRotationMatrix(phiB,thetaB,psiB ) ;
    
    % rotation matrix to strict body axis
    M1 = squeeze(rotM_roll(:,:,k)) * squeeze(rotM_YP(:,:,k)) ; 
    
    % pitch down by thetaB w.r.t body axis
    M2 = eulerRotationMatrix(0, -thetaB0rad, 0) ;
    M = M2 * M1 ;
    
    % M1' rotates a vector by (phi, theta, psi) in the lab frame, e.g.
    % M1' * [1;0;0] = data.AHat(k,:)'
    
    % M1  undoes the rotation of M1'. alternative description is that
    % M1 * data.AHat(k,:)' = [1;0;0]
    % give the coordinate of a lab-frame vector as described in the body frame
    
    % M2' * [1;0;0] = unit vector pitched down by thetab0
    % M2 * v = gived v pitched up by thetab0
    
    % first, represent the vector in the body-bound frame coordinates
    % then perform the rotation about the body roll axis
    % then represent back in whatever frame you want?
    
    % ---------------
    %% RIGHT WING
    % ---------------
    rotCR = M * (cr - cb) ; % rotated right wing center of mass
    rotSpanR = M * rightSpanHat ; % rotated right span
    rotChordR = M * rightChordHat ; % rotated right chord
    
    % get angles
    phiRdeg   = unwrap(atan2(rotSpanR(2),rotSpanR(1))) * RAD2DEG;
    thetaRdeg = asin(rotSpanR(3)) * RAD2DEG ;
    newEtaBody(k,1) = calcEta(rotSpanR, rotChordR,'right') ;
    etaRdeg = newEtaBody(k,1) ;
    % ---------------
    %% LEFT WING
    % ---------------
    rotCL = M * (cl - cb) ; % rotated left wing center of mass
    rotSpanL = M * leftSpanHat ; % rotated left span
    rotChordL = M * leftChordHat ; % rotated left chord
    
    % get angles
    phiLdeg   = unwrap(atan2(rotSpanL(2),rotSpanL(1))) * RAD2DEG;
    thetaLdeg = asin(rotSpanL(3)) * RAD2DEG ;
    newEtaBody(k,2) = calcEta(rotSpanL, rotChordL,'left') ;
    etaLdeg = newEtaBody(k,2) ;
    
    %------------------------------------------
    % store result
    anglesBodyFrame(k,:) = ...
        [0 0 phiRdeg thetaRdeg etaRdeg phiLdeg thetaLdeg etaLdeg] ;
    
    
end

%--------------------------------------------------------------------------
%% store angles in matrices

anglesBodyFrame(:,PHIR)   = + anglesBodyFrame(:,PHIR) ; % NOTE MINUS

anglesLabFrame  = unwrap(anglesLabFrame/RAD2DEG) * RAD2DEG ;
anglesBodyFrame = unwrap(anglesBodyFrame/RAD2DEG) * RAD2DEG ;

%--------------------------------------------------------------------------
%% plot results?
labels = {'\psi','beta','\phi_R','\theta_R', '\eta_R',...
    '\phi_L','\theta_L','\eta_L','\rho'} ;

plotReorder = [ PSI BETA PHIL PHIR ETAL ETAR THETAL THETAR ];
if (plotFlag)
    figure('name','Angles body frame') ;
    for s=1:8
        subplot(4,2,s) ;
        ind = plotReorder(s) ;
        plot(allT, anglesBodyFrame(:,ind),'ko-') ;
        xlabel('Time [ms]') ;
        title([labels{ind} ' body frame']) ;
        grid on ; box on ;
        set(gca,'xlim',[t(1) t(end)]);
    end
    figure('name','Angles lab frame') ;
    for s=1:8
        subplot(4,2,s) ;
        ind = plotReorder(s) ;
        plot(allT, anglesLabFrame(:,ind),'o-','color',[0 0.8 0]) ;
        xlabel('Time [ms]') ;
        title([labels{ind} ' lab frame']) ;
        grid on ; box on ;
        set(gca,'xlim',[allT(1) allT(end)]);
    end
    
end

end