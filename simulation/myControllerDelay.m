% -------------------------------------------------------------------------
% function to implement state-dependent delays for solving flight
% simulation dde
% -------------------------------------------------------------------------
function d = myControllerDelay(t, y, deltaT, omega)

% want the delay to correspond to a fixed time interval away from a
% fwd/back flip point (cycleTime = some fixed phase point in wingstroke)
cyclePhase = pi/2 ;
phaseCurr = wrapTo2Pi(rem(omega*t - cyclePhase, 2*pi)) ;  %omega*t - pi/2

cycleTime = t - phaseCurr./omega ;

fwdPhaseDiff = (3*pi)/2 - cyclePhase ;
fwdTimeDiff = fwdPhaseDiff/omega ;
d = cycleTime - (deltaT - fwdTimeDiff);

end



% cyclePhase1 = pi/2 ;
% cyclePhase2 = pi ;
% 
% % get phase of nearest current cycle point
% phaseCurr1 = wrapTo2Pi(rem(omega*t - cyclePhase1, 2*pi)) ;  %omega*t - pi/2
% phaseCurr2 = wrapTo2Pi(rem(omega*t - cyclePhase2, 2*pi)) ;  %omega*t - pi/2
% 
% % get times of nearest current cycle point
% cycleTime1 = t - phaseCurr1./omega ;
% cycleTime2 = t - phaseCurr2./omega ;
% 
% % get difference from fwd stroke to cycle times
% fwdPhaseDiff1 = (3*pi)/2 - cyclePhase1 ;
% fwdPhaseDiff2 = (3*pi)/2 - cyclePhase2 ;
% 
% fwdTimeDiff1 = fwdPhaseDiff1/omega ;
% fwdTimeDiff2 = fwdPhaseDiff2/omega ;
% 
% % d = [cycleTime1 - (deltaT - fwdTimeDiff1) ; ...
% %     cycleTime2 - (deltaT - fwdTimeDiff2)] ;
% d = [cycleTime1  ; ...
%     cycleTime2 ] ;