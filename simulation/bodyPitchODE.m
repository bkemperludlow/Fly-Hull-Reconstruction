function ds = bodyPitchODE(t,s) 
%% Using kinematics from Chang et al. (2014) to test different control models
%Notes: 
%   -theta in RADIANS!!!!
%   -state vector looks like: s = [y ydot z zdot theta thetadot]
%
%------------------------------------------------------------------------------
%% Initialize parameters

[~, omega, span, ~, body_mass, ~, MOI, g, airDensity,...
    S, r22_S, ~, R_sp, ~, hinge_vec, psi_m, psi_0, delta_psi, C,...
    phi_0, phi_m, K, ~, ~, ~, beta_0 ] = defineSimulationConstants() ;

% Define variables from state vector for convenience
%y = s(1) ;
%z = s(3) ; 
theta = s(5) ; %body pitch angle
thetadot = s(6) ; %rotational velocity

%% Define the Euler angles of the wing (set thetaR = thetaL = 0)

%Right wing:
phiR = phi_func(t) ;
%thetaL = 0 ;
psiR = psi_func(t) ; 

%Left wing:
phiL = phi_func(t) ;
%thetaL = 0 ;
psiL = psi_func(t) ;

%% Find wing motion using equations from Chang and Wang
% In the body frame
phiR_dot = (omega*(pi/180)*phi_m/asin(K))*cos(omega*t)./sqrt(1-K^2*(sin(omega*t)).^2) ; 
phiL_dot = (omega*(pi/180)*phi_m/asin(K))*cos(omega*t)./sqrt(1-K^2*(sin(omega*t)).^2) ; 

U_tR = phiR_dot*[span*cos((pi/180)*phiR),-span*sin((pi/180)*phiR), zeros(size(phiR))]' ;
U_tL = phiL_dot*[-span*cos((pi/180)*phiL),-span*sin((pi/180)*phiL), zeros(size(phiL))]' ;

%% Flip times
%{
%When the wing goes from forward to back stroke, not actual flip
temp = t*4*freq ; 
backFlipTimesInd = find(mod(temp,4) <= (1+1e-10) & mod(temp,4) >= (1-1e-10)) ;
fwdFlipTimesInd = find(mod(temp,4) <= (3+1e-10) & mod(temp,4) >= (3-1e-10)) ;
backFlipTimesR = t(backFlipTimesInd) ;
backFlipTimesL = t(backFlipTimesInd) ;
fwdFlipTimesR = t(fwdFlipTimesInd) ;
fwdFlipTimesL = t(fwdFlipTimesInd) ;
%}

%% Define angle of attack, so it can be used to calculate forces
% In the body frame (but not rotated by sp_angle)
rightSpanHat = [sin((pi/180)*phiR),cos((pi/180)*phiR), zeros(size(phiR))]' ;
leftSpanHat = [-sin((pi/180)*phiL),cos((pi/180)*phiL), zeros(size(phiL))]' ;
 
rightChordHat = [-cosd(psiR)*cosd(phiR); cosd(psiR)*sind(phiR) ; sind(psiR) ] ;
leftChordHat = [cosd(psiL)*cosd(phiL); cosd(psiL)*sind(phiL) ; sind(psiL)] ;

%rotates lab frame velocities into body coordinates
R = [ 1, 0, 0 ; 0, cos(pi/2-theta), (-1)*sin(pi/2-theta); ...
    0, sin(pi/2-theta), cos(pi/2-theta)] ;

%sum the velocities of the CoM and the wing (in body coordinates)
v_cm = [0 s(2) s(4)]' ;
thetadot_vec = [thetadot 0 0]' ; 
wingTipR = span*R_sp*rightSpanHat + hinge_vec ; 
wingTipL = span*R_sp*leftSpanHat + hinge_vec ; 

U_totR = R*v_cm + R_sp*U_tR + cross(thetadot_vec, wingTipR) ;
U_totL = R*v_cm + R_sp*U_tL + cross(thetadot_vec, wingTipL) ;
U_totR2 = dot(U_totR,U_totR) ;
U_totL2 = dot(U_totL,U_totL) ;

U_totR_hat = U_totR/sqrt(U_totR2) ; 
U_totL_hat = U_totL/sqrt(U_totL2) ; 

alphaR = (180/pi)*acos(dot(U_totR_hat,R_sp*rightChordHat,1)) ; 
alphaL = (180/pi)*acos(dot(U_totL_hat,R_sp*leftChordHat,1)) ; 

%% Calculate forces from the wings
%There are three coordinate systems that can be considered here.
%This function calculates the forces in the body frame by
%converting the center of mass velocity into body frame coordinates
%(just a rotation). Then it rotates the force vector back into lab
%coordinates

%drag and lift coefficients for left and right wings
C_L_R = .225 + 1.58*sin((pi/180)*(2.13*alphaR-7.2)) ;
C_D_R = 1.92 - 1.55*cos((pi/180)*(2.04*alphaR-9.82)) ;
C_L_L = .225 + 1.58*sin((pi/180)*(2.13*alphaL-7.2)) ;
C_D_L = 1.92 - 1.55*cos((pi/180)*(2.04*alphaL-9.82)) ;

%hat vectors for lift and drag (body coordinates)
LiftVec_R = cross(R_sp*rightSpanHat,U_totR_hat) ;
if LiftVec_R(3) < 0
    LiftVec_R = -LiftVec_R ;
end
DragVec_R = -U_totR_hat ;

LiftVec_L = cross(R_sp*leftSpanHat,U_totL_hat) ;
if LiftVec_L(3) < 0
    LiftVec_L = -LiftVec_L ;
end
DragVec_L = -U_totL_hat ;

%drag and lift for right wing
F_Lift_R = .5 * airDensity * S * r22_S * U_totR2 * C_L_R * LiftVec_R;
F_Drag_R = .5 * airDensity * S * r22_S * U_totR2 * C_D_R * DragVec_R ;

%drag and lift for left wing
F_Lift_L = .5 * airDensity * S * r22_S * U_totL2 * C_L_L * LiftVec_L ;
F_Drag_L = .5 * airDensity * S * r22_S * U_totL2 * C_D_L * DragVec_L ;

F_wing = F_Lift_R + F_Drag_R + F_Lift_L + F_Drag_L ; 
F_wing_lab = R'*F_wing ;

F_R_lab = R'*(F_Lift_R + F_Drag_R) ; 
F_L_lab = R'*(F_Lift_L + F_Drag_L) ; 

%% Find torque
%find lever arm vector in lab frame
CPR = .7*span*R_sp*rightSpanHat ; %centers of pressure are 70% along span
CPL = .7*span*R_sp*leftSpanHat ;
r_R = R'*(hinge_vec + CPR) ; %convert to lab frame
r_L = R'*(hinge_vec + CPL) ;

torqueR = cross(r_R, F_R_lab) ; 
torqueL = cross(r_L, F_L_lab) ; 
torqueTot = torqueR + torqueL ; 
% take only x component
torqueTot_x = dot(torqueTot,[1; 0; 0]) ;         

%% Set up differential equations
%state vector looks like s = [y, y_dot, z, z_dot, theta, theta_dot]
ds = zeros(6,1) ;

ds(1) = s(2) ;
ds(2) = (1/body_mass)*F_wing_lab(2) ; 

ds(3) = s(4) ; 
ds(4) = (1/body_mass)*F_wing_lab(3) - g ;
%disp(ds(4))

ds(5) = s(6) ;
ds(6) = (1/MOI)*torqueTot_x ;


%% Create data structure
%{
data = struct ;
%data.params.startTrackingTime = 0 ; %Need to change this
%data.params.endTrackingTime = N-1 ;
%data.params.fps = 1/dt ;
%data.params.N = N ;
%data.t = t ; 

data.params.span = span ;
data.params.chord = chord ;
data.params.freq = freq ;
data.params.psi_m = psi_m ;
data.params.delta_psi = delta_psi ;
data.params.psi_0 = psi_0 ;
data.params.pulseLengthMS = 0 ;
data.manualCorrRangeMS = [t1 t2]*1000 ;

data.rightWingTips = wingTipR/(50e-6) ; %convert to voxel units to match with 'data
data.leftWingTips = wingTipL/(50e-6) ;

data.fwdFlipTimesR = fwdFlipTimesR ;
data.backFlipTimesR = backFlipTimesR ;

data.fwdFlipTimesL = fwdFlipTimesL ;
data.backFlipTimesL = backFlipTimesL ;

data.rightSpanHats = rightSpanHats ;
data.leftSpanHats = leftSpanHats ;

data.rightChordHats = rightChordHats ; 
data.leftChordHats = leftChordHats ;

data.AHat = AHat ;
data.bodyCM = bodyCM ; %given in voxels so that it matches with 'data'

data.alphaR = alphaR ;
data.alphaL = alphaL ;

data.phiR = phiR ;
data.phiL = phiL ;

data.psiR = psiR ;
data.psiL = psiL ;

data.thetaR = thetaR ;
data.thetaL = thetaL ; 

data.U_tR = U_tR ;
data.U_tL = U_tL ;
%}
%% Misc functions
function phi = phi_func(t)
    phi = phi_0 + phi_m*asin(K*sin(omega*t))/asin(K) ;
end

function psi = psi_func(t)
    psi = psi_0 + psi_m*tanh(C*sin(omega*t + delta_psi))/tanh(C) ;
end

end


