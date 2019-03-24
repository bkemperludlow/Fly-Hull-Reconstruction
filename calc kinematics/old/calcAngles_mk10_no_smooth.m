function [anglesLabFrame, anglesBodyFrame, t, newEtaLab, newEtaBody, ...
    rho_fd_lambda, rho_fd_ci] = calcAngles_mk10_no_smooth(data, plotFlag) 
%
% this version uses a vector-based calculation of the wing pitch angles eta.
% this calculation is done in the function calcEta.m :
% calcEta(s, c, leftRight)
%
% wing angles in the body frame of reference are calculated differently.
% the body-F.O.R is defined as follows:
% let (xlab, ylab, zlab) be the lab frame of reference
% let (phib, thetab, psib) be the body yaw pitch and roll angles
% body axes are defined by rotating the lab axis using the following Euler
% rotation matrix:
% R(phi=phib, theta=NOMINAL THETA, psi=psib)
% NOMINAL THETA IS THE mean/representing TILT ANGLE of the stroke plane
% with respect to the body axis
%
% do not smooth to calculate angles. use the rho in the data structure
% which is already smoothed (given that it was obtained from the 3D hull
% rather than the pin).

% each matrix contains 8 columns for
% beta (body), psi (body), phiR, phiL, thetaR, thetaL, etaR, etaL

% subfunction used:
% calcBodyEulerAngles
% calcEta
disp('-----------------------------------------------------------------------') ;
disp('Using wing veins to calculate phi and theta for each wing');
disp('Using body-bound frame of reference with constant body pitch of theta_0')
disp('to calculate wing angles in the body frame') ;
disp('-----------------------------------------------------------------------') ;

if (~exist('plotFlag','var'))
    plotFlag = false ;
end


PSI    = 1 ; % body yaw
BETA   = 2 ; % body pitch
PHIR   = 3 ; 
THETAR = 4 ; 
ETAR   = 5 ; 
PHIL   = 6 ; 
THETAL = 7 ; 
ETAL   = 8 ; 
RHO    = 9 ;  %#ok<NASGU> % body roll


anglesLabFrame  = zeros(data.Nimages, 9);
anglesBodyFrame = zeros(data.Nimages, 8);

newEtaLab  = zeros(data.Nimages,2) ;
newEtaBody = zeros(data.Nimages,2) ;

% find roll angle - "true" Euler angle, i.e. rotation around the body axis
% can store the fd_lambda structure inside "data" such that the calc below
% will be done "only" once
fps = data.params.fps ;

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
C   = 180 / pi ;
deg2rad = pi / 180 ;

%tol = 1e-12 ;

%For now we assume the roll angle is zero, so we don't need lines 86-105. Need to update
rho0 = zeros(1,data.Nimages) ; 
%{
smoothRhoFlag = ~isfield(data,'rho_fd_lambda') ;

if (~smoothRhoFlag)
    smoothRhoFlag = isempty(data.rho_fd_lambda) ;
end

if (smoothRhoFlag)
    [~, smoothedRollAngle, rho_fd_lambda, rho_fd_ci ] = ...
        calcBodyEulerAngles (data, data.rhoTimes, plotFlag);
    
    rho0 = smoothedRollAngle ; % in degrees
else
    rho0 = zeros(1,data.Nimages) ;
    rho0(frameIndices) = eval_fd(t, data.rho_fd_lambda) ;
    
    rho_fd_lambda = [] ;
    rho_fd_ci     = [] ;
end
%}
% ----------------------------------------------------------
% calculate raw psi and raw beta for all frames (IN DEGREES)
% ----------------------------------------------------------
rawPsi  = zeros(data.Nimages,1) ;
rawBeta = zeros(data.Nimages,1) ;

for k=1:data.Nimages
    AHat = data.AHat(k,:) ;
    rawPsi(k)  = (180/pi)*atan2(AHat(2),AHat(1)); % body angle with respect to x axis IN DEGREES
    rawBeta(k) = (180/pi)*asin(AHat(3));  % body angle with respect to the horizon IN DEGREES
end
clear AHat;

psi0  = rawPsi ;
beta0 = rawBeta ;

thetaB0 = data.thetab0deg * deg2rad ;

for k=1:data.Nimages
    
    % ------------------
    % COPY RELEVANT DATA
    % ------------------
    cb = data.bodyCM(k,:)' ;
    cr = data.rightWingCM(k,:)' ;
    cl = data.leftWingCM(k,:)' ;
        
    rightChordHat = data.rightChordHats(k,:)' ;
    leftChordHat  = data.leftChordHats(k,:)' ;
    rightSpanHat  = data.rightSpanHats(k,:)' ;
    leftSpanHat   = data.leftSpanHats(k,:)' ;
    
%     rightVein = data.rightVeins(k,:)' ;    Changed by L&S. Need to update vein finding
%     leftVein  = data.leftVeins(k,:)' ; 
    rightVein = data.rightSpanHats(k,:)' ;
    leftVein  = data.leftSpanHats(k,:)' ; 
    
    newEtaLab(k,1) = calcEta(rightSpanHat, rightChordHat,'right') ;
    newEtaLab(k,2) = calcEta(leftSpanHat, leftChordHat,'left') ;
    
    etaRdeg = newEtaLab(k,1)  ;
    etaLdeg = newEtaLab(k,2)  ; 
    
    % Use the smoothed version of the roll angle rho
    bodyRollAngle = rho0(k) ; 
    
    % ---------------------------------------------------------------------
    % CALCULATE ANGLES IN THE LAB FRAME 
    % ---------------------------------------------------------------------
    
    % body  

    %psiRad  = psi0(k)*pi/180 ;    % body angle with respect to x axis
    %betaRad = beta0(k)*pi/180 ;   % body angle with respect to the horizon
    psiDeg  = psi0(k) ;
    betaDeg = beta0(k) ;
    
    % right wing
    phiRdeg   = atan2(rightVein(2),rightVein(1)) * C ;
    %phiRdeg   = myatan2pi(rightVein(2),rightVein(1)) * C ;    
    thetaRdeg = asin(rightVein(3)) * C;
    
    % left wing
    phiLdeg   = atan2(leftVein(2), leftVein(1)) * C;
    %phiLdeg   = myatan2pi(leftVein(2), leftVein(1)) * C;
    thetaLdeg = asin(leftVein(3)) * C;
    
    % store result
    anglesLabFrame(k,:) = ...
        [psiDeg betaDeg phiRdeg thetaRdeg etaRdeg phiLdeg...
        thetaLdeg etaLdeg bodyRollAngle] ;
      
    % -----------------------------------------------
    % CALCULATE ANGLES IN THE BODY FRAME OF REFERENCE
    % -----------------------------------------------
    
    % first rotation matrix bring (xlab, ylab, zlab) to 
    % (xbody ybody zbody) such that xbody is AHat
    
    %{
    phi   = -psiDeg*deg2rad ;
    % XXX
    theta = 0 ; % - betaDeg*deg2rad ; % - data.thetab0deg*deg2rad ; % constant theta
    psi   = -bodyRollAngle*deg2rad ;
    
    cph=cos(phi)    ; sph=sin(phi)   ;
    cth=cos(theta)  ; sth=sin(theta) ;
    cps=cos(psi)    ; sps=sin(psi)   ;
    
    lab_to_body1 = [cth*cph            cth*sph      (-sth) ; ...
        (sps*sth*cph-cps*sph) (sps*sth*sph+cps*cph) cth*sps ; ...
        (cps*sth*cph+sps*sph) (cps*sth*sph-sps*cph) cth*cps ] ;

    % second rotation matrix rotates (xbody ybody zbody) by -theta0 DOWN
    % about the ybody (pitch)
    phi = 0 ;
    % XXX
    theta = 0 ; % data.thetab0deg*deg2rad ; % constant theta (minus minus is plus)
    psi   = 0 ;
    
    cph=cos(phi)    ; sph=sin(phi)   ;
    cth=cos(theta)  ; sth=sin(theta) ;
    cps=cos(psi)    ; sps=sin(psi)   ;
    
    body1_to_body2 = [cth*cph            cth*sph           (-sth) ; ...
        (sps*sth*cph-cps*sph) (sps*sth*sph+cps*cph) cth*sps ; ...
        (cps*sth*cph+sps*sph) (cps*sth*sph-sps*cph) cth*cps ] ;

    M = ( body1_to_body2 * lab_to_body1)' ;
      %}
    
    phiB   = psiDeg*deg2rad ;
    thetaB = betaDeg*deg2rad ;
    psiB   = bodyRollAngle*deg2rad ;
    
    M1 = eulerRotationMatrix(phiB,thetaB,psiB ) ;        % to strict body axis
    
    %{
    %THIS PART SEEMS TO WORK
    %(this M2 is equivalent to calcAngles_mk9_no_smooth.m
    M2 = eulerRotationMatrix(0, -thetaB, 0) ;           % undo pitch only
        
    % the next M2 is a body-fixed F.O.R tilted by thetaB0
    %M2 = eulerRotationMatrix(0, -thetaB+thetaB0, 0) ;     % undo pitch only
    M3 =  eulerRotationMatrix(0, -(thetaB0-thetaB),0) ;
    
    M = M3 * M2 * M1 ;
    
    % M1*ahat = [1 0 0] ;
    % M1 * M2 * ahat =    body heading along x but pitching up by thetaB
    % M3 * M2 * M1 * ahat = body should pitch up by thetaB0
    %}
    
    M2 = eulerRotationMatrix(0, -thetaB0, 0) ; % pitch down by thetaB w.r.t body axis
    M = M2 * M1 ;
    
    % M1' rotates a vector by (phi, theta, psi) in the lab frame, e.g.
    % M1' * [1;0;0] = data.AHat(k,:)'
    
    % M1  undoes the rotation of M1'. alternative description is that
    % M1 * data.AHat(k,:)' = [1;0;0] 
    % give the coordinate of a lab-frame vector as described in the body frame

    
    % first, represent the vector in the body-bound frame coordinates
    % then perform the rotation about the body roll axis
    % then represent back in whatever frame you want?
    
    
    % RIGHT WING
    % ----------
    
    % right vein (-cb because everything is relative to body center of mass)
    rotHingeR = M * ( data.rightHinges(k,:)' - cb);
    rotTemp   = M * ( data.rightHinges(k,:)' + rightVein - cb);
    rotVeinR  = rotTemp - rotHingeR ;
    
    % right chord
    rotCR = M * (cr - cb) ;
    rotTemp   = M * ( cr + rightChordHat - cb) ;
    rotChordR = rotTemp - rotCR ;
    
    % right span
    rotTemp   = M * ( cr + rightSpanHat - cb) ;
    rotSpanR = rotTemp - rotCR ;
    
    % LEFT WING
    % ----------
    
    % left vein
    rotHingeL = M * ( data.leftHinges(k,:)' - cb);
    rotTemp   = M * (data.leftHinges(k,:)' + leftVein - cb) ;
    rotVeinL  = rotTemp - rotHingeL ;
    
    % left chord
    rotCL = M * (cl - cb) ;
    rotTemp   = M * ( cl + leftChordHat - cb) ;
    rotChordL = rotTemp - rotCL ;
    
    % left span
    rotTemp   = M * ( cl + leftSpanHat - cb) ;
    rotSpanL = rotTemp - rotCL ;
        

    % -----------------------------------------------------------------
    % CALCULATE PHI, THETA AND ETA FOR EACH WING IN THE BODY FRAME OF
    % REFERENCE, I.E. AFTER ROTATION.
    % -----------------------------------------------------------------
    
    % right wing
    % ----------
    
    phiRdeg   = atan2(rotVeinR(2),rotVeinR(1)) * C;
    thetaRdeg = asin(rotVeinR(3)) * C ;    
    newEtaBody(k,1) = calcEta(rotSpanR, rotChordR,'right') ;
    etaRdeg = newEtaBody(k,1) ;
    
    % left wing
    % ---------
    phiLdeg   = atan2(rotVeinL(2),rotVeinL(1)) * C;
    thetaLdeg = asin(rotVeinL(3)) * C ;
    newEtaBody(k,2) = calcEta(rotSpanL, rotChordL,'left') ;
    etaLdeg = newEtaBody(k,2) ;
    
    % store result
    anglesBodyFrame(k,:) = ...
        [0 0 phiRdeg thetaRdeg etaRdeg phiLdeg thetaLdeg etaLdeg] ;    
end


% ---------------------------------------
% ADJUST THE ANGLES TO THE CORRECT RANGES
% ---------------------------------------

%anglesBodyFrame(:,PHIR)   = + anglesBodyFrame(:,PHIR) ; % NOTE MINUS

anglesLabFrame  = unwrap(anglesLabFrame/C) * C ;
anglesBodyFrame = unwrap(anglesBodyFrame/C) * C ;



labels = {'\psi','beta','\phi_R','\theta_R', '\eta_R',...
    '\phi_L','\theta_L','\eta_L','\rho'} ;

plotReorder = [ PSI BETA PHIL PHIR ETAL ETAR THETAL THETAR ];

if (isfield(data,'ignoreFrames'))
    ignoreFrames = data.ignoreFrames ;
else
    ignoreFrames = [] ;
end

if (plotFlag)
    figure('name','Angles body frame') ;
    for s=1:8
        subplot(4,2,s) ;
        ind = plotReorder(s) ;
        plot(allT, anglesBodyFrame(:,ind),'ko-') ;
        hold on ;
        plot(allT, anglesBodyFrame(:,ind),'ko-') ;
        plot(allT(ignoreFrames), anglesBodyFrame(ignoreFrames,ind),'rx') ;
        xlabel('Time [sec]') ;
        title([labels{ind} ' body frame']) ;
        grid on ; box on ;        
        set(gca,'xlim',[t(1) t(end)]);
    end    
    figure('name','Angles lab frame') ;
    for s=1:8
        subplot(4,2,s) ;
        ind = plotReorder(s) ;
        hold on ;
        plot(allT, anglesLabFrame(:,ind),'o-','color',[0 0.8 0]) ;
        plot(allT(ignoreFrames), anglesLabFrame(ignoreFrames,ind),'rx') ;
        xlabel('Time [sec]') ;
        title([labels{ind} ' lab frame']) ;
        grid on ; box on ;
        set(gca,'xlim',[allT(1) allT(end)]);
    end
    hold off;
end

end