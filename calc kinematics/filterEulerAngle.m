%--------------------------------------------------------------------------
% function to smooth the body/wing Euler angles from a fly video using a
% low-pass butterworth filter. 
%
% NB: writing this and similar functions to make the smoothing process
% consistent across code base
%--------------------------------------------------------------------------
function angle_filt = filterEulerAngle(angle_raw, filter_lvl)
% -------------------------------------------------------
%% filter out jumps at endpoints that lead to transients
% alternative: use padarray to replicate signal on both sides
angle_raw_hampel = hampel(angle_raw) ; 
% -------------------------------
%% apply low pass filter
d1 = designfilt('lowpassiir','FilterOrder',3,'SampleRate',8000, ...
    'HalfPowerFrequency',filter_lvl,'DesignMethod','butter'); %hpf = 100
angle_filt = filtfilt(d1,angle_raw_hampel) ;

end