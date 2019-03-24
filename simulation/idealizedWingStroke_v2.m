function data = idealizedWingStroke_v2(phi_front,phi_back)
%Want to generate a fake set of wing data to test my torque calculations,
%and later my control model fit. 
%
%--------------------------------------------------------------------------
%% Define some of the aspects of the wing flight, e.g. amplitude, frequency.
%Later, I'll want these to be functions, so I should try to keep it general

freq = 250 ; %Hz. This doesn't seem to shift between wings
omega = 2*pi*freq ;

psi_m = 53 ; %deg. This the amplitude for wing pitch
psi_0 = 90 ; 
delta_psi = 72.4*(pi/180) ; %rad. Wang uses -72.4 degrees
C = 2.4 ;

phi_m = (phi_back - phi_front)/2 ; %deg (this is the stroke amplitude). was 155
phi_0 = (phi_front + phi_back)/2 ;
K = .7 ;

theta_m = 25 ; %guesses [deg]
theta_0 = 0 ;
delta_theta = 90*(pi/180) ; %70.6*(pi/180) ; %from Berman

span = .0025 ; %meters
chord = .0007 ; %meters
%body_mass = 1.1e-6 ; %kg
%body_length = .0024 ; %meters
%body_width = .0012 ; %meters
%MOI = .506e-12 ; %N*m^2

r_hinge = [0 .00018 .0006] ; %m

t1 = 0 ;
t2 = .02 ;
dt = 1e-5 ;
t = t1:dt:t2 ;
N = length(t) ;

%% Define the Euler angles of the wing

%Right wing:
phiR =  phi_0 + phi_m*asin(K*sin(omega*t))/asin(K) ; 
thetaR = theta_0 + theta_m*cos(2*omega*t + delta_theta) ;
%alphaR = alpha_0 + beta*cos(2*pi*freq*t + phaseDiff) ;
psiR = psi_0 + psi_m*tanh(C*sin(omega*t + delta_psi))/tanh(C) ; %This needs to be changed

%Left wing:
phiL = phi_0 + phi_m*asin(K*sin(omega*t))/asin(K) ; 
thetaL = theta_0 + theta_m*cos(2*omega*t + delta_theta) ;
%alphaL = alpha_0 + beta*cos(2*pi*freq*t + phaseDiff);
psiL = psi_0 + psi_m*tanh(C*sin(omega*t + delta_psi))/tanh(C) ; %This needs to be changed

%% Find wing motion using equations from Wang et al.

rightSpanHats = [sin((pi/180)*phiR).*cos((pi/180)*thetaR); cos((pi/180)*phiR).*cos((pi/180)*thetaR);...
    sin((pi/180)*thetaR)]' ;
leftSpanHats = [-sin((pi/180)*phiL).*cos((pi/180)*thetaL); cos((pi/180)*phiL).*cos((pi/180)*thetaR);...
    sin((pi/180)*thetaL)]' ;
 
%rightChordHats = [-cosd(psiR).*cosd(phiR); cosd(psiR).*sind(phiR) ; sind(psiR) ]' ;
%leftChordHats = [cosd(psiL).*cosd(phiL); cosd(psiL).*sind(phiL) ; sind(psiL)]' ;
%only works for xy stroke plane!
rightChordHats = [-(cos((pi/180)*psiR).*cos((pi/180)*phiR) + sin((pi/180)*psiR).*sin((pi/180)*phiR).*sin((pi/180)*thetaR)) ;...
    sin((pi/180)*phiR).*cos((pi/180)*psiR) - sin((pi/180)*psiR).*cos((pi/180)*phiR).*sin((pi/180)*thetaR) ;...
    cos((pi/180)*thetaR).*sin((pi/180)*psiR) ]' ;
leftChordHats = [cos((pi/180)*psiL).*cos((pi/180)*phiL) + sin((pi/180)*psiL).*sin((pi/180)*phiL).*sin((pi/180)*thetaL) ;...
    sin((pi/180)*phiL).*cos((pi/180)*psiL) - sin((pi/180)*psiL).*cos((pi/180)*phiL).*sin((pi/180)*thetaL) ;...
    cos((pi/180)*thetaL).*sin((pi/180)*psiL) ]' ;

wingTipR = span*rightSpanHats ;
wingTipL = span*leftSpanHats ;

AHat = repmat([0, 1/sqrt(2), 1/sqrt(2)],N,1) ;
bodyCM = zeros(N,3) - repmat(r_hinge,N,1) ; %given in voxels so that it matches with 'data'


%% Wing flip times

temp = t*4*freq ; 
backFlipTimesInd = find(mod(temp,4) <= (1+1e-10) & mod(temp,4) >= (1-1e-10)) ;
fwdFlipTimesInd = find(mod(temp,4) <= (3+1e-10) & mod(temp,4) >= (3-1e-10)) ;
backFlipTimesR = t(backFlipTimesInd) ;
backFlipTimesL = t(backFlipTimesInd) ;
fwdFlipTimesR = t(fwdFlipTimesInd) ;
fwdFlipTimesL = t(fwdFlipTimesInd) ;


%% Need to make a structure like 'data' so that I can feed this into quasisteadyTorqueSam

data = struct ;
data.params.startTrackingTime = 0 ; %Need to change this
data.params.endTrackingTime = N-1 ;
data.params.fps = 1/dt ;
data.params.N = N ;
data.t = t ; 
data.params.span = span ;
data.params.chord = chord ;
data.params.freq = freq ;

data.params.psi_0 = psi_0 ;
data.params.delta_psi = delta_psi ;
data.params.psi_m = psi_m ;

data.params.phi_0 = phi_0 ;
data.params.phi_m = phi_m ;

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
data.bodyCM = bodyCM/(50e-6) ; %given in voxels so that it matches with 'data'

%data.alphaR = alphaR ;
%data.alphaL = alphaL ;

data.phiR = phiR ;
data.phiL = phiL ;

data.psiR = psiR ;
data.psiL = psiL ;

data.thetaR = thetaR ;
data.thetaL = thetaL ; 

%data.thetaDot = 0 ;
%{
if phi_front <= 20 
    data.thetaDot = -3000*(pi/180) ; %rad/s
else
    data.thetaDot = +3000*(pi/180) ;
end
%}
%data.U_tR = U_tR ;
%data.U_tL = U_tL ;

%save idealizedWingstrokeData_phasediff data
end
   