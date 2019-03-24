function [anglesLabFrame, anglesBodyFrame, t, newEtaLab, newEtaBody, sp_rho, smoothed_rho, rho_t, rho_samp] = ...
    calcAngles_quick_and_dirty_mk3(data, rollErr, plotFlag)
%
% this version uses a vector-based calculation of the wing pitch angles eta.
% this calculation is done in the function calcEta.m :
% calcEta(s, c, leftRight)
%
% do not smooth to calculate angles. use the rho in the data structure
% which is already smoothed (given that it was obtained from the 3D hull
% rather than the pin).

% each matrix contains 8 columns for
% beta (body), psi (body), phiR, phiL, thetaR, thetaL, etaR, etaL

% subfunction used:
% calcBodyEulerAngles
% calcEta


%{
data2 = data ;
data2.rightVeins = rightVeins ;
data2.leftVeins  = leftVeins ;
data2.rightHinges = rightHinges ;
data2.leftHinges  = leftHinges ;
%}
disp('NOT Using wing veins to calculate phi and theta for each wing');

if (~exist('plotFlag','var'))
    plotFlag = false ;
end


PSI    = 1 ;
BETA   = 2 ;
PHIR   = 3 ;
THETAR = 4 ;
ETAR   = 5 ;
PHIL   = 6 ;
THETAL = 7 ;
ETAL   = 8 ;
RHO    = 9 ; %#ok<NASGU>


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
tol = 1e-12 ;

% set roll angle to zero. quick and dirty.
%rho0 = zeros(1, data.Nimages) ;


% ----------------------------------------------------------
% calculate raw psi and raw beta for all frames (IN DEGREES)
% ----------------------------------------------------------
rawPsi  = zeros(data.Nimages,1) ; % yaw
rawBeta = zeros(data.Nimages,1) ; % pitch

for k=1:data.Nimages
    AHat = data.AHat(k,:) ;
    rawPsi(k)  = (180/pi)*atan2(AHat(2),AHat(1)); % body angle with respect to x axis IN DEGREES
    rawBeta(k) = (180/pi)*asin(AHat(3));  % body angle with respect to the horizon IN DEGREES
end


% calc roll angles at data.rhoTimes
if (~isfield(data,'rhoTimes')) || (isempty(data.rhoTimes))
    %[rhoTimes, rollVectors] = estimateRollVectors(data);
    disp('data.rhoTimes is not defined.') ;
    [rhoTimes, rollVectors] = estimateRollVector(data) ;
else
    rhoTimes    = data.rhoTimes ;
    rollVectors = data.rollVectors ;
end


if (rhoTimes(1)~=1)
    rollVectors(1,:) = rollVectors(rhoTimes(1),:) ;
    %rollVectors = [ rollVectors(1,:) ; rollVectors] ; % when using estimateRollVectos
    rhoTimes = [1 rhoTimes] ;
end
if (rhoTimes(end)~=data.Nimages)
    rollVectors(end,:) = rollVectors(rhoTimes(end),:) ;
    %rollVectors(end+1,:) = rollVectors(end,:) ; % when using estimateRollVectos
    rhoTimes = [rhoTimes data.Nimages] ; % when using estimateRollVectos
end

C = 180 / pi ;
rho_samp = zeros(length(rhoTimes),1) ;
rho_samp2 = rho_samp ;

rho_t    = rhoTimes + data.params.firstTrackableFrame - 1 ;
rho_t    = rho_t / data.params.fps  ;
for j=1:length(rhoTimes)
    it = rhoTimes(j) ;
    xb = data.AHat(it, :) ;
    thetab = rawBeta(it) ;  % body pitch in degrees
    phib   = rawPsi(it) ;   % body yaw in degrees.
    
    yb = rollVectors(it,:) ; %  goes from R-->L wing hinges
    
    % unpitch matrix
    unpitch = eulerRotationMatrix(0, -thetab*pi/180, 0) ;
    
    % unyaw matrix
    unyaw = eulerRotationMatrix(-phib*pi/180, 0, 0) ;
    
    rotyb = unpitch' * unyaw' * yb' ;
    
    rho_samp(j) = acos( rotyb(2) ) * C * sign(rotyb(3)) ;
    
    %{
    % previous way - it is equivalent to the one above
    %yb = rollVectors(j,:) ; %  goes from R-->L wing hinges % when using estimateRollVectos
    zb = cross(xb,yb) ;
    
    %yaw   = atan2(xb(2), xb(1)) * C ;
    %pitch = asin(xb(3)) * C ;
    
    xb_proj = [xb(1) xb(2) 0 ] ;
    xb_proj = xb_proj / norm(xb_proj) ;
    v = cross([0 0 1], xb_proj) ;
    
    rho_samp2(j) = acos( dot(v, yb)) * C * sign(yb(3));
    %}
end

%{
figure ;
plot(rho_samp,'bo-') ;
hold on ;
plot(rho_samp2,'r+-');
keyboard ;
%}

ESTERR = rollErr ; % 4 ;
%dt     = mean(diff(rhoTimes)) / data.params.fps ; % because rho is sampled every ~half stroke
%currN  = length(rhoTimes) ;

tol = ESTERR ;

[sp_rho, ~] =  mySplineSmooth(rho_t, rho_samp, tol) ;
rho0 = fnval(sp_rho, allT) ;

smoothed_rho = rho0 ;


    

%{
figure ;
plot(rho_t*1000, rho_samp,'ko','markerfacecolor','g') ;
hold on ;
plot(allT*1000, rho0,'k-') ;
xlabel('Time [ms]') ;
ylabel('Roll angle [deg]') ;
%}
%keyboard

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

psi0  = rawPsi ;
beta0 = rawBeta ;

if (~isfield(data,'thetaB0'))
    thetaB0rad = 45 * pi / 180 ;
    disp('data.thetaB0 is not defined. taking 45 degrees as default value') ;
else
    thetaB0rad = data.thetaB0 * pi / 180 ;
end

deg2rad = pi / 180 ;

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
    
    %rightVein = data.rightVeins(k,:)' ;
    %leftVein  = data.leftVeins(k,:)' ;        
    rightVein = rightSpanHat ; % data.rightVeins(k,:) ;
    leftVein  = leftSpanHat ; % data.leftVeins(k,:) ;
    
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
    
    M2 = eulerRotationMatrix(0, -thetaB0rad, 0) ; % pitch down by thetaB w.r.t body axis
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
    
    
    % RIGHT WING
    % ----------
    
    % right vein (-cb because everything is relative to body center of mass)
    %rotHingeR = M * ( data.rightHinges(k,:)' - cb);
    %rotTemp   = M * ( data.rightHinges(k,:)' + rightVein - cb);
    %rotVeinR  = rotTemp - 0 ;
    rotVeinR = M * rightVein ;
    
    
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
    %rotHingeL = M * ( data.leftHinges(k,:)' - cb);
    %rotTemp   = M * (data.leftHinges(k,:)' + leftVein - cb) ;
    %rotVeinL  = rotTemp - rotHingeL ;
    rotVeinL = M * leftVein ;
    
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


%{

for k=1:data.Nimages
    
    % ------------------
    % COPY RELEVANT DATA
    % ------------------
    cb = data.bodyCM(k,:) ;
    cr = data.rightWingCM(k,:) ;
    cl = data.leftWingCM(k,:) ;
    
    %AHat = data.AHat(k,:) ;
    
    rightChordHat = data.rightChordHats(k,:) ;
    leftChordHat  = data.leftChordHats(k,:) ;
    rightSpanHat  = data.rightSpanHats(k,:) ;
    leftSpanHat   = data.leftSpanHats(k,:) ;
    
    rightVein = rightSpanHat ; % data.rightVeins(k,:) ;
    leftVein  = leftSpanHat ; % data.leftVeins(k,:) ;
    
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
    
    % First rotation - fix body roll angle by rotating (-bodyRollAngle)
    % degrees about the BODY axis AHat
    
    AHat = data.AHat(k,:) ;
    
    % axis of rotation in real space is the line from p1 to p2 ;
    % NOTE that here roll is defined along the BODY axis (Euler angle)
    p1 = cb ;
    p2 = cb + AHat ;
    angle = -bodyRollAngle ;
    
    % using p1=[0 0 0] and p2=AHat should work as well,
    % because only p2-p1 matter
    
    % no point in rotating AHat, since it stays the same - the rotation
    % axis itself is AHat
    % rotAHat = RotatePoint( data.AHat(k,:), p1, p2, angle) ;
    
    % rotate right wing c.m.
    rotCR     = cb + RotatePoint (cr-cb, p1, p2, angle) ;
    
    % rotate right chord vector
    chordTemp = cb + RotatePoint (cr+rightChordHat-cb, p1, p2, angle) ;
    rotChordR = chordTemp - rotCR ;
    
    % rotate right span vector
    spanTemp  = cb + RotatePoint (cr+rightSpanHat-cb, p1, p2, angle) ;
    rotSpanR  = spanTemp - rotCR ;
    
    % rotate right vein
    % taking cr and rotCR in the following two lines gives the same results
    % as taking rightWingHinge and the rotated rightWingHinge
    veinTemp  = cb + RotatePoint (cr+rightVein-cb, p1, p2, angle) ;
    rotVeinR  = veinTemp - rotCR ;
    
    
    % rotate left wing c.m.
    rotCL     = cb + RotatePoint (cl-cb, p1, p2, angle) ;
    
    % rotate left chord vector
    chordTemp = cb + RotatePoint (cl+leftChordHat-cb, p1, p2, angle) ;
    rotChordL = chordTemp - rotCL ;
    
    % rotate left span vector
    spanTemp  = cb + RotatePoint (cl+leftSpanHat-cb, p1, p2, angle) ;
    rotSpanL  = spanTemp - rotCL ;
    
    % rotate left vein
    veinTemp  = cb + RotatePoint (cl+leftVein-cb, p1, p2, angle) ;
    rotVeinL  = veinTemp - rotCL ;
    
    clear spanTemp chordTemp veinTemp
    
    % SECOND ROTATION
    % ---------------
    % rotate the whole fly such that its heading is along
    % the x axis. This simplifies the calculation of angles in the body
    % frame of reference.
    
    % next rotation axis points in parallel the z axis
    p1 = cb ;
    p2 = cb + [ 0 0 1 ] ;
    
    angle = -psiDeg ;
    
    % rotate body axis
    
    rotAHat = RotatePoint( AHat, p1, p2, angle) ;
    
    psi2Deg  = atan2(rotAHat(2),rotAHat(1)) * C ;  % new psi should be zero
    beta2Deg = asin(rotAHat(3)) * C;               % new beta should be the same as old beta
    
    
    % verify body axis rotation
    if ( abs(psi2Deg)>tol )
        disp ('new psi is not zero - check this...') ;
        keyboard ;
    end
    
    if (abs(beta2Deg-betaDeg)>tol)
        disp('new beta is different than original beta - check this...') ;
        keyboard ;
    end
    
    % rotate right wing c.m.
    rotCR2     = cb + RotatePoint (rotCR-cb, p1, p2, angle) ;
    
    % rotate right chord vector
    chordTemp  = cb + RotatePoint (rotCR+rotChordR-cb, p1, p2, angle) ;
    rotChordR2 = chordTemp - rotCR2 ;
    
    % rotate right span vector
    spanTemp  = cb + RotatePoint (rotCR+rotSpanR-cb, p1, p2, angle) ;
    rotSpanR2  = spanTemp - rotCR2 ;
    
    % rotate right vein
    veinTemp  = cb + RotatePoint (rotCR+rotVeinR-cb, p1, p2, angle) ;
    rotVeinR2  = veinTemp - rotCR2 ;
    
    % rotate left wing c.m.
    rotCL2     = cb + RotatePoint (rotCL-cb, p1, p2, angle) ;
    
    % rotate left chord vector
    chordTemp = cb + RotatePoint (rotCL+rotChordL-cb, p1, p2, angle) ;
    rotChordL2 = chordTemp - rotCL2 ;
    
    % rotate left span vector
    spanTemp  = cb + RotatePoint (rotCL+rotSpanL-cb, p1, p2, angle) ;
    rotSpanL2  = spanTemp - rotCL2 ;
    
    % rotate left vein
    veinTemp  = cb + RotatePoint (rotCL+rotVeinL-cb, p1, p2, angle) ;
    rotVeinL2  = veinTemp - rotCL2 ;
    
    clear spanTemp chordTemp tempVein
    
    % -----------------------------------------------------------------
    % RECALCULATE PHI, THETA AND ETA FOR EACH WING IN THE BODY FRAME OF
    % REFERENCE, I.E. AFTER ROTATION.
    % -----------------------------------------------------------------
    
    % right wing
    % ----------
    
    phiRdeg   = atan2(rotVeinR2(2),rotVeinR2(1)) * C;
    thetaRdeg = asin(rotVeinR2(3)) * C ;
    newEtaBody(k,1) = calcEta(rotSpanR2, rotChordR2,'right') ;
    etaRdeg = newEtaBody(k,1) ;
    
    % left wing
    % ---------
    phiLdeg   = atan2(rotVeinL2(2),rotVeinL2(1)) * C;
    thetaLdeg = asin(rotVeinL2(3)) * C ;
    newEtaBody(k,2) = calcEta(rotSpanL2, rotChordL2,'left') ;
    etaLdeg = newEtaBody(k,2) ;
    
    % store result
    anglesBodyFrame(k,:) = ...
        [psi2Deg beta2Deg phiRdeg thetaRdeg etaRdeg phiLdeg thetaLdeg etaLdeg] ;
end
%}

% ---------------------------------------
% ADJUST THE ANGLES TO THE CORRECT RANGES
% ---------------------------------------

anglesBodyFrame(:,PHIR)   = + anglesBodyFrame(:,PHIR) ; % NOTE MINUS

anglesLabFrame  = unwrap(anglesLabFrame/C) * C ;
anglesBodyFrame = unwrap(anglesBodyFrame/C) * C ;

%{
% adjust body frame RIGHT phi
ind = find(anglesBodyFrame(:,PHIR)>90) ;
anglesBodyFrame(ind,PHIR) = anglesBodyFrame(ind,PHIR) - 360 ;
anglesBodyFrame(:,PHIR)   = - anglesBodyFrame(:,PHIR) ; % NOTE MINUS
% adjust body frame LEFT phi
ind = anglesBodyFrame(:,PHIL) < -90 ;
anglesBodyFrame(ind,PHIL) = 360 + anglesBodyFrame(ind,PHIL) ;

% adjust lab frame RIGHT phi
%{
ind = find(anglesLabFrame(:,PHIR)>90) ;
anglesLabFrame(ind,PHIR) = anglesLabFrame(ind,PHIR) - 360 ;
ind = find(anglesLabFrame(:,PHIR)<-90) ;
anglesLabFrame(ind,PHIR) = 360 + anglesLabFrame(ind,PHIR) ;
anglesLabFrame(:,PHIR) = - anglesLabFrame(:,PHIR) ;
%}

% crap
%anglesLabFrame(:,PHIR) = mod(anglesLabFrame(:,PHIR), 360) ;
%anglesLabFrame(:,PHIL) = mod(anglesLabFrame(:,PHIL), 360) ;

%{
% adjust body frame RIGHT eta
ind = find(anglesBodyFrame(:,ETAR)<-90) ;
anglesBodyFrame(ind,ETAR) = 360 + anglesBodyFrame(ind,ETAR) ;
% adjust body frame LEFT eta
ind = find(anglesBodyFrame(:,ETAL)<-90) ;
anglesBodyFrame(ind,ETAL) = 360 + anglesBodyFrame(ind,ETAL) ;
%}

% adjust body frame RIGHT eta
ind = find(anglesBodyFrame(:,ETAR)>300) ;
anglesBodyFrame(ind,ETAR) = 360 - anglesBodyFrame(ind,ETAR) ;

% adjust body frame RIGHT eta
ind = find(anglesBodyFrame(:,ETAL)>300) ;
anglesBodyFrame(ind,ETAL) = 360 - anglesBodyFrame(ind,ETAL) ;

%}



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