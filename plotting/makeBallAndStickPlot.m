% make wing motion diagram 

defineConstantsScript
wingType = 'L' ;

span = 2.5 ; %mm
chord = 0.7/2 ; %mm

t = (data.params.startTrackingTime : data.params.endTrackingTime ) / 8000 ;
if strcmp(wingType,'R')
    fwdFlipTimes = data.fwdFlipTimesR ;
    backFlipTimes = data.backFlipTimesR ;
    phi = -data.anglesBodyFrame(:,PHIR) ;
    theta = data.anglesBodyFrame(:,THETAR) ;
    psi = data.anglesBodyFrame(:,ETAR) ;
elseif strcmp(wingType,'L')
    fwdFlipTimes = data.fwdFlipTimesL ;
    backFlipTimes = data.backFlipTimesL ;
    phi = data.anglesBodyFrame(:,PHIL) ;
    theta = data.anglesBodyFrame(:,THETAL) ;
    psi = data.anglesBodyFrame(:,ETAL) ;
else
    keyboard ;
end

phiEstErr = 1 ;
thetaEstErr = 1 ; 
psiEstErr = 2.5 ;

[sp_phi, ~, ~] = mySplineSmooth(t, phi, phiEstErr) ;
[sp_theta, ~, ~] = mySplineSmooth(t, theta, phiEstErr) ;
[sp_psi, ~, ~] = mySplineSmooth(t, psi, phiEstErr) ;

phi_smoothed = fnval(sp_phi,t) ;
theta_smoothed = fnval(sp_theta, t) ; 
psi_smoothed = fnval(sp_psi,t) ; 


dot_x = span*cos((pi/180)*phi_smoothed) ;
dot_y = span*sin((pi/180)*theta_smoothed) ;

stick_x = chord*cos((pi/180)*psi_smoothed) ; 
stick_y = chord*sin((pi/180)*psi_smoothed) ; 

startFlipInd = 11 ; 

i1 = find(t == fwdFlipTimes(startFlipInd)) ;
i2 = find(t == backFlipTimes(startFlipInd+1)) ;
i3 = find(t == fwdFlipTimes(startFlipInd+1)) ;
indSpacing = 1 ;

figure ;
subplot(2,1,1)
hold on
plot(dot_x(i1:indSpacing:i2), dot_y(i1:indSpacing:i2), 'ko','MarkerFaceColor','k')
for k = i1:indSpacing:i2
    plot([dot_x(k) (dot_x(k) - stick_x(k))], [dot_y(k) (dot_y(k) - stick_y(k))],'k-') 
end
axis equal

subplot(2,1,2)
hold on
plot(dot_x(i2:indSpacing:i3), dot_y(i2:indSpacing:i3), 'ko','MarkerFaceColor','k')
for k = i2:indSpacing:i3
    plot([dot_x(k) (dot_x(k) - stick_x(k))], [dot_y(k) (dot_y(k) - stick_y(k))],'k-') 
end
axis equal


