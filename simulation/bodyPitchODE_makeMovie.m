function [] = bodyPitchODE_makeMovie(t,y,z,theta, savePath)
%% Visualize the fly dynamics in 3D (basically performs the same functions again)
%
%
%
%---------------------------------------------------------------------------------
[~, omega, span, chord, ~, body_length, ~, ~, ~,...
    ~, ~, ~, R_sp, ~, hinge_vec, psi_m, psi_0, delta_psi, C,...
    phi_0, phi_m, K, ~, ~, ~, ~, ~, ...
    ~, ~, ~] = defineSimulationConstants() ;

N = length(t) ;

hfly = figure ;
ax = [-.002 .002 -.003 .003 -.003 .002] ;

cd(savePath) ; 
imgCounter = 0 ; 

for i = 1:N
    imgCounter = imgCounter + 1 ; 
    figure(hfly) ;
    clf ;
    hold on ;
    
    phiR = phi_0 + phi_m*asin(K*sin(omega*t(i)))/asin(K) ;
    psiR = psi_0 + psi_m*tanh(C*sin(omega*t(i) + delta_psi))/tanh(C) ;
    
    phiL = phi_0 + phi_m*asin(K*sin(omega*t(i)))/asin(K) ;
    psiL = psi_0 + psi_m*tanh(C*sin(omega*t(i) + delta_psi))/tanh(C) ;
    
    
    wingTipR = [span*sin((pi/180)*phiR),span*cos((pi/180)*phiR), zeros(size(phiR))]' ;
    wingTipL = [-span*sin((pi/180)*phiL),span*cos((pi/180)*phiL), zeros(size(phiL))]' ;
    
    rightChord = chord*[-cosd(psiR)*cosd(phiR); cosd(psiR)*sind(phiR) ; sind(psiR) ] ;
    leftChord = chord*[cosd(psiL)*cosd(phiL); cosd(psiL)*sind(phiL) ; sind(psiL)] ;
    
    %rotates lab frame velocities into body coordinates
    R = [ 1, 0, 0 ; 0, cos(pi/2-theta(i)), (-1)*sin(pi/2-theta(i)); ...
        0, sin(pi/2-theta(i)), cos(pi/2-theta(i))] ;
    
    centerOfMass = [0, y(i), z(i)]' ;
    
    bodyAxis = R'*[0 0 body_length]' ;
    bodyBack = R'*[0 0 -body_length/2]' + centerOfMass ;
    hingeLab = R'*hinge_vec ;
    
    myPlotVector(bodyBack,bodyAxis,'k',4) ;
    myPlotVector(centerOfMass,hingeLab,'k',3) ;
    
    %Plot wing spans
    wingBase = hingeLab + centerOfMass ;
    rightWingTip_lab = R'*R_sp*wingTipR ;
    leftWingTip_lab = R'*R_sp*wingTipL ;
    
    myPlotVector(wingBase,rightWingTip_lab,'r',3.5) ;
    myPlotVector(wingBase,leftWingTip_lab,'b',3.5) ;
    
    %Plot wing chords
    chordBaseR = wingBase + .5*rightWingTip_lab ;
    chordBaseL = wingBase + .5*leftWingTip_lab ;
    rightChord_lab = R'*R_sp*rightChord ;
    leftChord_lab = R'*R_sp*leftChord ;
    
    myPlotVector(chordBaseR,rightChord_lab,[0 .5 0],2.5) ;
    myPlotVector(chordBaseL,leftChord_lab,[.5 0 .5],2.5) ;
    
    tstr = num2str(t(i),4) ;
    
    title(['t = ' tstr])
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    set(gca,'fontsize',12)
    
    hold off ;
    axis equal ;
    grid on ;
    box on ;
    
    az = mod(111 + i,360)  ;
    el = 26 ;
    view([az,el]) ;
    
    axis(ax)
    %axis tight
    print(gcf,['pitch_simulation_' num2str(imgCounter)],'-dpng','-r0') ;
    disp(t(i))
    pause(.5) ; 
end
