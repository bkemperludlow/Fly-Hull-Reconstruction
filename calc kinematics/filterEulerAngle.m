%--------------------------------------------------------------------------
% function to smooth the body/wing Euler angles from a fly video using a
% low-pass butterworth filter. 
%
% NB: writing this and similar functions to make the smoothing process
% consistent across code base
%--------------------------------------------------------------------------
function angle_filt = filterEulerAngle(angle_raw, filter_lvl)

d1 = designfilt('lowpassiir','FilterOrder',3,'SampleRate',8000, ...
    'HalfPowerFrequency',filter_lvl,'DesignMethod','butter'); %hpf = 100
angle_filt = filtfilt(d1,angle_raw) ;

end