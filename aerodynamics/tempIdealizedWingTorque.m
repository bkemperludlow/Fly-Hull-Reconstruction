
phi_front = [0:55]' ; 
runSize = length(phi_front) ;
phi_back = 175*ones(runSize,1) ;
timeAvgTorque = zeros(runSize,1) ;
sumTorque = zeros(runSize,1) ;
maxTorque = zeros(runSize,1) ;

MOI = .506e-12 ;

for i = 1:runSize
    clear('data')
    data = idealizedWingStroke_v2(phi_front(i),phi_back(i)) ;
    run quasiSteadyTorqueSam.m
    
    [~,j1] = min(abs(t - data.backFlipTimesR(2))) ;
    [~,j2] = min(abs(t - data.backFlipTimesR(2+1))) ;
    sumTorque(i) = sum(data.pitchTorque(j1:j2)) ;
    
    beatLength = data.backFlipTimesR(2+1) - data.backFlipTimesR(2) ;
    dt = 1/data.params.fps ;
    timeAvgTorque(i) = sumTorque(i)*dt/beatLength ;
    
    [maxNorm,j3] = max(abs(data.pitchTorque(j1:j2))) ;
    maxTorque = maxNorm*sign(data.pitchTorque(j1+j3)) ;
    
end

%{
save phi_front phi_front
save phi_back phi_back
save timeAvgTorque timeAvgTorque
save sumTorque sumTorque
%}
