dataPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Controller Analysis\' ;
cd(dataPath)

load accelerationAndPhiStruct
%load controllerAnalysisStruct

Nmovies = length(accelerationAndPhiStruct) ;
controlInd = find([accelerationAndPhiStruct(:).flyType] == 2) ;
experimentalInd = find([accelerationAndPhiStruct(:).flyType] == 1) ;
deltaPhiFrontCell = cell(Nmovies,1) ; 
pitchVelCell = cell(Nmovies,1) ; 
pitchVelMax = zeros(Nmovies,1) ; 
deltaPitchMax = zeros(Nmovies,1) ; 
deltaPitchCell = cell(Nmovies,1) ; 
deltaPhiFrontMax = zeros(Nmovies,1) ; 

for i = 1:Nmovies 
    
    %if i == 14 || i == 15
    %    continue ; 
    %end
    fwdFlipTimes = accelerationAndPhiStruct(i).fwdFlipTimes ;
    backFlipTimes = accelerationAndPhiStruct(i).backFlipTimes ;
    phiFront = accelerationAndPhiStruct(i).phiFront ;
    bodyPitch = accelerationAndPhiStruct(i).bodyPitch ;
    t = accelerationAndPhiStruct(i).t ;
    pitchVel = accelerationAndPhiStruct(i).pitchVelocity ;
    flyType = accelerationAndPhiStruct(i).flyType ;
    peakInd = accelerationAndPhiStruct(i).peakInd ;
    
    pitchVel_avg = zeros((length(backFlipTimes) - 1),1) ;
    deltaPitch_avg = zeros((length(backFlipTimes) - 1),1) ;
    fwdFlipTimes_new = zeros((length(backFlipTimes) - 1),1) ;
    deltaPhiFront = zeros((length(backFlipTimes) - 1),1) ;
    
    bodyPitchPrePert = mean(bodyPitch(t <= 0 & t > -0.01)) ;
    phiFrontPrePert = mean(phiFront(fwdFlipTimes > -0.01 & fwdFlipTimes < 0)) ;
    
    deltaPitch = bodyPitch - bodyPitchPrePert ;
    %maxPert = abs(deltaPitch(peakInd)) ; 
    %maxCorrection = abs(accelerationAndPhiStruct(i).deltaPhiFront) ;
    
    for j = 1:(length(backFlipTimes) - 1)
        [~,i1] = min(abs(t - backFlipTimes(j))) ;
        [~,i2] = min(abs(t - backFlipTimes(j+1))) ;
        pitchVel_avg(j) = mean(pitchVel(i1:i2)) ;
        deltaPitch_avg(j) = mean(deltaPitch(i1:i2))  ;
        
        fwdFlipTimeInd = find(fwdFlipTimes < backFlipTimes(j+1) & fwdFlipTimes > backFlipTimes(j), 1, 'first') ; 
        fwdFlipTimes_new(j) = fwdFlipTimes(fwdFlipTimeInd) ;
        deltaPhiFront(j) = phiFront(fwdFlipTimeInd) - phiFrontPrePert ;
        
        if abs(deltaPhiFront(j)) > 62
            keyboard ; 
        end
    end
    
    [~, pitchVelPeakLocs] = findpeaks(abs(pitchVel(t > 0))) ; 
    [~,deltaPhiFrontPeakLocs] = findpeaks(abs(deltaPhiFront(fwdFlipTimes_new > 0))) ; 
    t0ind = find(t == 0) ; 
    t0FlipInd = find(fwdFlipTimes_new < 0,1, 'last') ;
    if isempty(t0FlipInd)
        t0FlipInd = 0 ;
    end
    pitchVelMax(i) = pitchVel(pitchVelPeakLocs(1)+t0ind) ;
    deltaPitchMax(i) = deltaPitch(peakInd) ; 
    try
        deltaPhiFrontMax(i) = deltaPhiFront(deltaPhiFrontPeakLocs(1)+t0FlipInd) ; 
    catch
        keyboard ;
    end
    
    if (0)
        figure ; 
        subplot(2,1,1)
        hold on
        plot(t,deltaPitch,'k-')
    end
    
    pitchVelCell{i} = pitchVel_avg ;
    deltaPhiFrontCell{i} = deltaPhiFront ;
    deltaPitchCell{i} = deltaPitch_avg  ;
    
end

h_vel = figure ; 
hold on
for k = 1:Nmovies
    deltaPhiFront_temp = deltaPhiFrontCell{k} ;
    pitchVel_temp = pitchVelCell{k} ;
    
    if ~isempty(find(controlInd == k,1))
        color_temp = [0 .7 0] ;
    elseif ~isempty(find(experimentalInd == k,1))
        color_temp = [0.7 0 0] ;
    else
        continue ; 
    end
    
    plot(deltaPhiFront_temp(2:end),pitchVel_temp(1:end-1),'ko','MarkerFaceColor',color_temp)
end
xlabel('\Delta \Phi_{front}')
ylabel('$\dot{\theta}$','interpreter','latex')

controlDeltaPhiFront = [] ;
controlDeltaPitch = [] ;
experimentalDeltaPhiFront = [] ;
experimentalDeltaPitch = [] ; 
h_disp = figure ; 
hold on
for k = 1:Nmovies
    deltaPhiFront_temp = deltaPhiFrontCell{k} ;
    deltaPitch_temp = deltaPitchCell{k} ;
    
    if ~isempty(find(controlInd == k,1))
        color_temp = [0 .7 0] ;
        controlDeltaPhiFront = [controlDeltaPhiFront ; deltaPhiFront_temp] ;
        controlDeltaPitch = [controlDeltaPitch ; deltaPitch_temp] ;
        %continue ; 
    elseif ~isempty(find(experimentalInd == k,1))
        color_temp = [0.7 0 0] ;
        experimentalDeltaPhiFront = [experimentalDeltaPhiFront ; deltaPhiFront_temp] ;
        experimentalDeltaPitch = [experimentalDeltaPitch ; deltaPitch_temp] ;
        
    else
        continue ; 
    end
    
    plot(deltaPhiFront_temp(2:end),deltaPitch_temp(1:end-1),'ko','MarkerFaceColor',color_temp)
end
xlabel('\Delta \Phi_{front}')
ylabel('\Delta \theta')

h_both = figure ;
hold on
for k = 1:Nmovies
    deltaPhiFront_temp = deltaPhiFrontCell{k} ;
    deltaPitch_temp = deltaPitchCell{k} ;
    pitchVel_temp = pitchVelCell{k} ;
    
    if ~isempty(find(controlInd == k,1))
        color_temp = [0 .7 0] ;
    elseif ~isempty(find(experimentalInd == k,1))
        color_temp = [0.7 0 0] ;
    else
        continue ; 
    end
    
    plot3(deltaPitch_temp(1:end-1),pitchVel_temp(1:end-1),deltaPhiFront_temp(2:end),'ko','MarkerFaceColor',color_temp)
end
xlabel('\Delta \theta')
ylabel('$\dot{\theta}$','interpreter','latex')
zlabel('\Delta \Phi_{front}')
axis tight
box on 
grid on

figure ; 
hold on
plot(deltaPitchMax(controlInd), pitchVelMax(controlInd), 'ksq','MarkerFaceColor',[0 .7 0])
plot(deltaPitchMax(experimentalInd), pitchVelMax(experimentalInd), 'ksq','MarkerFaceColor',[.7 0 0])
xlabel('\Delta \theta_{max}')
ylabel('$\dot{\theta}_{max}$','interpreter','latex')

figure ; 
hold on
plot(deltaPitchMax(controlInd), deltaPhiFrontMax(controlInd), 'ksq','MarkerFaceColor',[0 .7 0])
plot(deltaPitchMax(experimentalInd), deltaPhiFrontMax(experimentalInd), 'ksq','MarkerFaceColor',[.7 0 0])
xlabel('\Delta \theta_{max}')
ylabel('$\Delta\phi_{max}$','interpreter','latex')

figure ; 
hold on
plot(pitchVelMax(controlInd), deltaPhiFrontMax(controlInd), 'ksq','MarkerFaceColor',[0 .7 0])
plot(pitchVelMax(experimentalInd), deltaPhiFrontMax(experimentalInd), 'ksq','MarkerFaceColor',[.7 0 0])
xlabel('$\dot{\theta}_{max}$','interpreter','latex')
ylabel('$\Delta\phi_{max}$','interpreter','latex')

%{
h_corr_vel = figure ;
hold on
for k = 1:Nmovies
    deltaPhiFront_temp = deltaPhiFrontCell{k} ;
    %deltaPitch_temp = deltaPitchCell{k} ;
    pitchVel_temp = pitchVelCell{k} ;
    
    if ~isempty(find(controlInd == k,1))
        color_temp = [0 .7 0] ;
    elseif ~isempty(find(experimentalInd == k,1))
        color_temp = [0.7 0 0] ;
    else
        continue ; 
    end
    
    [r, lags] = xcorr(pitchVel_temp,deltaPhiFront_temp,4,'coeff') ;
    plot(lags,r,'Color',color_temp)
end
title('Velocity Corr')

h_corr_pitch = figure ;
hold on
for k = 1:Nmovies
    deltaPhiFront_temp = deltaPhiFrontCell{k} ;
    deltaPitch_temp = deltaPitchCell{k} ;
    %pitchVel_temp = pitchVelCell{k} ;
    
    if ~isempty(find(controlInd == k,1))
        color_temp = [0 .7 0] ;
    elseif ~isempty(find(experimentalInd == k,1))
        color_temp = [0.7 0 0] ;
    else
        continue ;
    end
    
    [r, lags] = xcorr(deltaPitch_temp,deltaPhiFront_temp,4,'coeff') ;
    plot(lags,r,'Color',color_temp)
end
title('Angle Corr')
%}