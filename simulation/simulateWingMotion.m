%Runs a lot of 'idealizedWingStroke' calculations for different wing
%kinematics, and reports the torque

phi_front = 0:1:65 ;
phi_back = 175 ;
MOI = 7.17e-13 ; %kg*m^2

runSize = length(phi_front) ; 
pitchTorqueAvg = zeros(runSize,1) ;

%hpitchtorque = figure ;
%hold on
%patchColor = [1 1 1]*.8 ;

for i = 1:runSize
    data = idealizedWingStroke(phi_front(i),phi_back) ; 
    run quasiSteadyTorqueSam
    
    pulseLength = data.params.pulseLengthMS ;
    t = data.t ;
    N = length(t) ;
    tms = 1000*t ;
    dt = 1/data.params.fps ;
    
    fwdFlipTimesR = data.fwdFlipTimesR ;
    backFlipTimesR = data.backFlipTimesR ;
    fwdFlipTimesL = data.fwdFlipTimesL ;
    backFlipTimesL = data.backFlipTimesL ;
    
    %torqueR = data.torqueR ;
    %torqueL = data.torqueL ;
    torques = data.torques ;
    axisHats = data.AHat ;
    
    bodyY = cross(axisHats,repmat([0 0 1],N,1)) ;
    bodyY = bodyY ./ repmat(myNorm(bodyY),1,3) ;
    
    %pitchTorqueR = dot(torqueR,bodyY,2) ;
    %pitchTorqueL = dot(torqueL,bodyY,2) ;
    pitchTorque = dot(torques,bodyY,2) ;
    
    %temp1 = abs(t - backFlipTimesR(1)) ; %need to fix this
    %[~,i1] = min(temp1) ; 
    %temp2 = abs(t - fwdFlipTimesR(1)) ; %need to fix this
    %[~,i2] = min(temp2) ;

    %pitchTorqueAvg(i) = dt*sum(pitchTorque(i1:i2))/(t(i2)-t(i1)) ; 
    [maxPitchTorque,maxInd] = max(pitchTorque) ;
    [minPitchTorque,minInd] = min(pitchTorque) ;
    pitchTorqueAvg(i) = maxPitchTorque + minPitchTorque ;
    
    %currColor = [i/65 0 0] ;
    %if mod(i,5) == 0 
    %    plot(1000*t,pitchTorque,'Color',currColor)
    %end
end

%ylim = get(gca,'ylim');
%plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
%set(gca,'ylim',ylim) ;


hsimulate = figure ; 
hold on
set(gcf,'paperPositionMode','auto')
set(gca,'fontsize',14)

plot(phi_front, (1e-5)*(180/pi)*(pitchTorqueAvg-pitchTorqueAvg(21))/MOI, 'k.')
ylabel('Pitch Acceleration [\times10^5 deg/s^2]');
xlabel('Front Stroke Angle [deg]');
title('Toy Model Calculation')

