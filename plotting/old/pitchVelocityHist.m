%Checks splines in a general way and then plots histogram of maximum pitch
%velocity (counts are per individual)

splineCheckFlag = true ;

Nmovies = length(noPertBodyAngles) ;
nbins = 10 ;
%--------------------------------------------------------------------------------
%Plot time series for all events
if splineCheckFlag
    hsplinecheck = figure ;
    hold on
    for i =1:Nmovies
        
        if noPertBodyAngles(i).TNT
            plotColor = [.7 0 0 ] ;
        elseif ~noPertBodyAngles(i).TNT
            plotColor = [ 0 .7 0] ;
        else
            plotColor = [0 0 0] ;
            disp('No experiment number found')
        end
        
        smoothed_pitch = fnval(noPertBodyAngles(i).PitchSpline, noPertBodyAngles(i).Time) ;
        tms = 1000*noPertBodyAngles(i).Time ;
        plot(tms, smoothed_pitch,'Color',plotColor)
    end
    xlabel('time [ms]')
    ylabel('Pitch angle [deg]')
    axis tight
end
%--------------------------------------------------------------------------------

%Get maximum pitch velocity (just doing the easiest thing first, could use
%refining)
maxPitchVelTNT = nan(Nmovies,1) ;
maxPitchVelIMPTNT = nan(Nmovies,1) ;

stdPitchVelTNT = nan(Nmovies,1) ;
stdPitchVelIMPTNT = nan(Nmovies,1) ;

for j =1:Nmovies 
    sp_pitch = noPertBodyAngles(j).PitchSpline ;
    t = noPertBodyAngles(j).Time ;
    pitchVel = fnval( fnder(sp_pitch,1), t) ;
    
    [maxAbsPitchVel, Ind] = max( abs(pitchVel)) ;
        
    if noPertBodyAngles(j).TNT
        maxPitchVelTNT(j) = pitchVel(Ind) ;
        stdPitchVelTNT(j) = std(pitchVel) ;
    elseif ~noPertBodyAngles(j).TNT
        maxPitchVelIMPTNT(j) = pitchVel(Ind) ;
        stdPitchVelIMPTNT(j) = std(pitchVel) ;
    else
        disp('ERROR: Do not know if control or experimental group')
    end
end

[NMaxCtrl, xMaxCtrl ] = hist(maxPitchVelIMPTNT,nbins) ;
[NMaxExp, xMaxExp ] = hist(maxPitchVelTNT,nbins) ;
xlim1 = [min(xMaxExp) max(xMaxExp)] ;


[NStdCtrl, xStdCtrl ] = hist(stdPitchVelIMPTNT,nbins) ;
[NStdExp, xStdExp ] = hist(stdPitchVelTNT,nbins) ;
xlim2 = [min(xStdExp) max(xStdExp)] ;
%--------------------------------------------------------------------------------------
%plot histogram of maxima
hMaxVelDist = figure ;

s11 = subplot(2,1,1) ;
b11 = bar(xMaxCtrl, NMaxCtrl,'edgecolor','k','facecolor',[0 .7 0]) ;
set(gca,'xlim',xlim1)
xlabel('Max Pitch Velocity [deg/s]')
ylabel('Counts (# Flies)')
%title('Maximum Pitch Velocity (magnitude) Distribution')

s12 = subplot(2,1,2) ;
b12 = bar(xMaxExp, NMaxExp,'edgecolor','k','facecolor',[.7 0 0]) ;
set(gca,'xlim',xlim1)
xlabel('Max Pitch Velocity [deg/s]')
ylabel('Counts (# Flies)')

%--------------------------------------------------------------------------------------
%plot histogram of pitch velocity standard deviation (is this even relevant?)


hStdVelDist = figure ;

s21 = subplot(2,1,1) ;
b21 = bar(xStdCtrl, NStdCtrl,'edgecolor','k','facecolor',[0 .7 0]) ;
set(gca,'xlim',xlim2)
xlabel('Pitch Velocity Standard Deviation [deg/s]')
ylabel('Counts (# Flies)')
%title('Pitch Velocity Standard Deviation Distribution')

s22 = subplot(2,1,2) ;
b22 = bar(xStdExp, NStdExp,'edgecolor','k','facecolor',[.7 0 0]) ;
set(gca,'xlim',xlim2)
xlabel('Pitch Velocity Standard Deviation [deg/s]')
ylabel('Counts (# Flies)')

