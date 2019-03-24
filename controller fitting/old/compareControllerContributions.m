function [I_norm, P_norm] = ...
    compareControllerContributions(data, controller_fit_struct)

%get relevant data
defineConstantsScript
bodyPitch = data.anglesLabFrame(:,BETA) ;
K_i = controller_fit_struct.K_i ;
K_p = controller_fit_struct.K_p ;
deltaT = controller_fit_struct.deltaT ;

t = (data.params.startTrackingTime : data.params.endTrackingTime ) / 8000 ;
dt = 1/8000 ; 

if isfield(data,'manualCorrRangeMS')
    t_start = max([data.manualCorrRangeMS(1) -10])/1000 ; 
    t_end = min([data.manualCorrRangeMS(2) 30])/1000 ;
    t_range = t_start : dt : t_end  ; %seconds!
else 
    t_range = -0.01 : dt : 0.03 ; %seconds!
end

i1 = find(t == t_range(1)) ; 
i2 = find(t == t_range(end)) ;

%smooth pitch
d1 = designfilt('lowpassiir','FilterOrder',8,'SampleRate',8000, ...
    'HalfPowerFrequency',100,'DesignMethod','butter');
pitch_filt = filtfilt(d1,bodyPitch) ;
c_pitch = fit(t',pitch_filt,'linearinterp');

%pitch_smoothed = c_pitch(t) ;
pitchDisp = c_pitch(t-deltaT)- c_pitch(0) ;
pitchVel = differentiate(c_pitch,t-deltaT) ;

pitchDispL2 = sqrt(trapz(t_range, pitchDisp(i1:i2).^2)) ; 
pitchVelL2 = sqrt(trapz(t_range, pitchVel(i1:i2).^2)) ; 
%proportionalTerm = K_p*pitchVel ; 
%integralTerm = K_i*(c_pitch(t-deltaT)- c_pitch(0)) ; 

I_norm = K_i*pitchDispL2 ;
P_norm = K_p*pitchVelL2 ;


end