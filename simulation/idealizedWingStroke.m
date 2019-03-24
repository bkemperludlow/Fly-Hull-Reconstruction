function data = idealizedWingStroke(phi_front,phi_back,alpha_0,phaseDiff)
%Want to generate a fake set of wing data to test my torque calculations,
%and later my control model fit. So this should try to stay general.
%
%General conventions:
%   phi = wing stroke angle 
%   theta = wing stroke plane deviation angle 
%   psi = wing pitch angle 
%   alpha = angle of attack (defined as the angle between the wing chord
%       and the velocity of far-field flow)
%   beta = body pitch angle 
%   phaseDiff = phase difference between wingstroke and wing flip
%
%To simulate the wing flip we'll use the model from Wang et al. (2004):
%
%   alpha = alpha_0 + beta*sin(2*pi*freq + phaseDiff) 
%
%--------------------------------------------------------------------------
%% Define some of the aspects of the wing flight, e.g. amplitude, frequency.
%Later, I'll want these to be functions, so I should try to keep it general

freq = 200 ; %Hz. This doesn't seem to shift between wings
beta = 45 ; %deg
if (~exist('alpha_0','var'))
    alpha_0 = 90 ; %default is 90
end
if (~exist('phaseDiff','var'))
    phaseDiff = 0 ; %Wang et al (2004) use pi/4, 0, and -pi/4 to model advanced, symmetric, and delayed wing rotation
end

Amp = phi_back - phi_front ; %deg (this is the stroke amplitude). was 155
AmpR = Amp ;
AmpL = Amp ;

midStroke = .5*(phi_front + phi_back) ; %deg. was 97.5
midStrokeR = midStroke ;
midStrokeL = midStroke ;

t1 = 0 ;
t2 = .03 ;
dt = 1e-5 ;
t = t1:dt:t2 ;
N = length(t) ;

%% Define the Euler angles of the wing

%Right wing:
phiR = (AmpR/2)*sin(2*pi*freq*t) + midStrokeR ;
thetaR = 0 ;
alphaR = alpha_0 + beta*cos(2*pi*freq*t + phaseDiff) ;
psiR = alphaR ; %This needs to be changed

%Left wing:
phiL = (AmpL/2)*sin(2*pi*freq*t) + midStrokeL ;
thetaL = 0 ;
alphaL = alpha_0 + beta*cos(2*pi*freq*t + phaseDiff);
psiL = alphaL ; %This needs to be changed

%% Find wing motion using equations from Wang et al.

span = .0025 ; %in meters
chord = .0007 ; %in meters

wingTipR = [span*sin((pi/180)*phiR)',span*cos((pi/180)*phiR)', zeros(size(phiR))'] ;
wingTipL = [-span*sin((pi/180)*phiL)',span*cos((pi/180)*phiL)', zeros(size(phiL))'] ;

U_tR = (AmpR*pi^2*freq/180)*[span*cosd(phiR).*cos(2*pi*freq*t); -span*sind(phiR).*cos(2*pi*freq*t); zeros(size(phiR))]' ;
U_tL = (AmpL*pi^2*freq/180)*[-span*cosd(phiL).*cos(2*pi*freq*t); -span*sind(phiL).*cos(2*pi*freq*t); zeros(size(phiL))]' ;

rightSpanHats = [sind(phiR)',cosd(phiR)', zeros(size(phiR))'] ;
leftSpanHats = [-sind(phiL)',cosd(phiL)', zeros(size(phiL))'] ;
 
rightChordHats = [-cosd(alphaR).*cosd(phiR); cosd(alphaR).*sind(phiR) ; sind(alphaR) ]' ;
leftChordHats = [cosd(alphaL).*cosd(phiL); cosd(alphaL).*sind(phiL) ; sind(alphaL)]' ;

AHat = repmat([0, 1/sqrt(2), 1/sqrt(2)],N,1) ;
bodyCM = zeros(N,3) ; %given in voxels so that it matches with 'data'
r_hinge = [0 -.0002 0] ; % meters SIGN HERE???
bodyCM = bodyCM - (1/50e-6)*repmat(r_hinge,N,1) ; 

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
data.params.beta = beta ;
data.params.phaseDiff = phaseDiff ;
data.params.alpha_0 = alpha_0 ;
data.params.pulseLengthMS = 0 ;
data.manualCorrRangeMS = [t1 t2]*1000 ;

data.rightWingTips = wingTipR/(50e-6) ; %convert to voxel units to match with 'data
data.leftWingTips = wingTipL/(50e-6) ;

data.fwdFlipTimesR = fwdFlipTimesR ;
data.backFlipTimesR = backFlipTimesR ;

data.fwdFlipTimesL = fwdFlipTimesL ;
data.backFlipTimesL = backFlipTimesL ;

data.backFlipTimesInd = backFlipTimesInd ;
data.fwdFlipTimesInd = fwdFlipTimesInd ;

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

%save idealizedWingstrokeData_phasediff data