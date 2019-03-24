%Given the torques on the fly (which I calculate), I should be able to
%qualitatively capture the dynamics

%Definition of terms:
%   I = moment of inertia (about pitch axis)
%   alpha = angular acceleration
%   omega = angular velocity 
%   beta = pitch angle
%
%
%---------------------------------------------------------------------
%Preliminaries for finding data, etc.
defineConstantsScript
exprNum = 7 ;
movNum  = 65 ;

if (movNum<10)
    zstr = '00' ;
elseif (movNum<100)
    zstr = '0' ;
else
    zstr = '';
end

datapath = ['F:\luca\Analysis\pitch down\Expr_' num2str(exprNum) '_mov_' zstr num2str(movNum) '\' ] ;
cd(datapath)

datafilename = [ datapath ...
    'Expr' num2str(exprNum) 'mov' zstr num2str(movNum) '_Data_manually_corrected.mat' ] ; %

load(datafilename) ;

%Define constants
dt = 1/8000 ; %1/fps
I = .506e-12 ; %from Cheng et al. (2009). units are N*m*s^2 = kg*m^2

%Data that needs to be loaded in
t = data.t ;
torques = data.torques ; 
axisHats = data.AHat ;
N = length(t) ;

%Preallocate arrays
alpha = zeros(N,1) ;
omega = zeros(N,1) ;
beta = zeros(N,1) ;

%define inital and final times over which to calculate torque
t1 = .00625 ; %in seconds. must be multiple of 1/8000
t2 = .02 ;

N1 = find(t == t1) ;
N2 = find(t == t2) ;

%Find magnitude of torque in the pitch direction
bodyY = cross(axisHats,repmat([0 0 1],N,1)) ;
bodyY = bodyY ./ repmat(myNorm(bodyY),1,3) ;

pitchTorque = dot(torques,bodyY,2) ;

%Find real pitch angle
realPitch = (180/pi)*asin(axisHats(:,3)) ;

%{
pitchEstErr = .5 ;
[sp_pitch, pitch_smooth, ~] = mySplineSmooth(t,realPitch,pitchEstErr) ;

%check spline fit

hsmooth = figure ;
plot(t,realPitch,'b.')
hold on
plot(t, fnval(sp_pitch,t),'r')
%}

%Define inital values
alpha(N1) = (180/pi)*(1/I)*pitchTorque(N1) ;
omega(N1) = ((180/pi)*asin(axisHats(N1+1,3))-(180/pi)*asin(axisHats(N1-1,3)))/dt ; %fnval( fnder(sp_pitch,1), t1) ; 
beta(N1) = (180/pi)*asin(axisHats(N1,3)) ;

%The dynamics
for i = (N1+1):N2
    alpha(i) = (180/pi)*(1/I)*pitchTorque(i) ;
    omega(i) = omega(i-1) + alpha(i-1)*dt ;
    beta(i) = beta(i-1) + omega(i-1)*dt + .5*alpha(i-1)*(dt)^2 ;
end

%create figure to examine differences

hpitchcomp = figure ;
plot(t(N1:N2),beta(N1:N2),'b','lineWidth',2)
hold on
plot(t(N1:N2), realPitch(N1:N2), 'r', 'lineWidth',2)
axis tight

