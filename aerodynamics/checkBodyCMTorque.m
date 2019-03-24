phiFront = linspace(0, 40, 41) ;
phiBack = linspace(140,180,41) ;

[phiF, phiB] = meshgrid(phiFront,phiBack) ;
sumTorque = zeros(size(phiF)) ;

N = length(phiFront)*length(phiBack) ;

for i = 1:N
    data = idealizedWingStroke(phiF(i),phiB(i)) ;
    run quasiSteadyTorqueSam
    i1 = data.backFlipTimesInd(3) ;
    i2 = data.backFlipTimesInd(4) ;
    sumTorque(i) = sum(data.torques(i1:i2,1)) ;
    disp(i)
end

figure ; surf(phiF,phiB,abs(sumTorque))
colorbar

avgCos = besselj(0,(pi/180)*(phiB-phiF)/2).*cos((pi/180)*(phiB+phiF)/2) ;
diff = abs(p_2/R + avgCos) ;
figure ; surf(phiF,phiB,diffCheck) ;