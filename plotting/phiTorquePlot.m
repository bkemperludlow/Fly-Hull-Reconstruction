function [h, timeAvgTorque, thetaErr] = phiTorquePlot(data) 
% NEED TO HAVE DATA LOADED AND BE IN THE PROPER DIRECTORY FOR
% THIS VERSION OF THE CODE. WILL FIX LATER

titleStr = 'Quasi-Steady Torque Calculation' ;

useOldPlot = 0 ;

try
    pulseLength = data.params.pulseLengthMS ;
catch
    pulseLength = 0 ;
end
try
    t = data.t ;
catch
    t1 = data.params.startTrackingTime ;
    t2 = data.params.endTrackingTime ;
    dt = 1/data.params.fps ;
    t = t1:dt:t2 ;
end
N = length(t) ;
tms = 1000*t ;

fwdFlipTimesR = data.fwdFlipTimesR ;
backFlipTimesR = data.backFlipTimesR ;
fwdFlipTimesL = data.fwdFlipTimesL ;
backFlipTimesL = data.backFlipTimesL ;


try
    pitchTorque = data.pitchTorque ;
catch
    %torqueR = data.torqueR ;
    %torqueL = data.torqueL ;
    torques = data.torques ;
    axisHats = data.AHat ;

    bodyY = cross(axisHats,repmat([0 0 1],N,1)) ;
    bodyY = bodyY ./ repmat(myNorm(bodyY),1,3) ;

    %pitchTorqueR = dot(torqueR,bodyY,2) ;
    %pitchTorqueL = dot(torqueL,bodyY,2) ;
    pitchTorque = dot(torques,bodyY,2) ;
end
patchColor = [1 1 1 ] * 0.8 ; 
faceAlpha = 1 ;
fps = data.params.fps ; 
dt = 1/fps ; 


tsfvec = [0 pulseLength pulseLength 0 0] ;
if isfield(data, 'manualCorrRangeMS')
    xlim = data.manualCorrRangeMS ;
elseif isfield(data, 'correctionTime')
    xlim = data.correctionTime ;
else 
    xlim = [tms(1) tms(end)] ;
end
i1 = find(tms == xlim(1)) ; %i1 and i2 give the frame limits of the manual correction
i2 = find(tms >= (xlim(2)-(1e-9)),1,'first') ;
%ylim = [-3 3]*1e-9 ;

h = figure ;
    %set(gcf, 'PaperPositionMode', 'manual');
    %set(gcf, 'PaperUnits', 'inches');
    %set(gcf, 'PaperPosition', [2 1 9 3]);
    hold on;
    
    htemp = plot(tms(i1:i2), pitchTorque(i1:i2)*1e8,'Color',[.5 0 0], 'LineWidth', 2) ;
    ylim = get(gca,'ylim') ;
    delete(htemp) ; 
    
    set(gca,'fontsize',10) ;
    avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 238 170]/255,'facealpha',faceAlpha) ;
    plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
    
    hcalc = plot(tms, pitchTorque*1e8,'Color',[0 .5 0], 'LineWidth', 2) ;
    harea = area(tms, pitchTorque*1e8,'LineStyle','none','FaceColor',[0 .5 0]) ;
    plot(tms, zeros(size(tms)),'k-','LineWidth',.5)
    %legend([hf,hcalc],{'Pert.','\tau_{total}'},'location','northwest') ;
    %grid on ; box on ;
    ylabel('Pitch Torque [\times10^{-8} N*m]') ;
    set(gca,'ylim',ylim ) ;
    set(gca,'xlim',xlim ) ;
    xlabel('Time [ms]')
    title(titleStr) ;

if nargout > 1
    Ntemp = length(backFlipTimesR)-1 ;
    torqueErr = zeros(Ntemp,1) ;
    timeAvgTorque = zeros(Ntemp,1) ;
    thetaErr = zeros(Ntemp,1) ;
    thetaDotErr = zeros(Ntemp,1) ;
    %thetaErr2 = zeros(Ntemp,1) ;
    MOI = .506e-12 ;
    for i = 1:Ntemp
        [~,i1] = min(abs(t - backFlipTimesR(i))) ;
        [~,i2] = min(abs(t - backFlipTimesR(i+1))) ;
        torqueErr(i) = sum(pitchTorque(i1:i2)) ;
        beatLength = backFlipTimesR(i+1) - backFlipTimesR(i) ;
        timeAvgTorque(i) = torqueErr(i)*dt/beatLength ;
        thetaChange = (180/pi)*cumtrapz(t(i1:i2),cumtrapz(t(i1:i2),pitchTorque(i1:i2)))/MOI ;
        thetaDotChange = (180/pi)*cumtrapz(t(i1:i2),pitchTorque(i1:i2))/MOI ; 
        thetaErr(i) = thetaChange(end) ;
        thetaDotErr(i) = thetaDotChange(end) ;
    end
    disp('Time-averaged torque = ')
    disp(timeAvgTorque)
    disp('Cumulative angle change = ')
    thetaErr
    disp('Cumulative anglular velocity change = ')
    thetaDotErr
end
%{
hPhiTorque = figure('position',[140 50 1100 900]) ;

s1 = subplot(3,1,1); %Total torques
    hold on;
    
    htemp = plot(tms(i1:i2), pitchTorque(i1:i2),'Color',[.5 0 0], 'LineWidth', 2) ;
    ylim = get(gca,'ylim') ;
    delete(htemp) ; 
    
    set(gca,'fontsize',14) ;
    avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
    
    plot(tms, pitchTorque,'Color',[0 .5 0], 'LineWidth', 2) ;
    legend({'Pert.','\tau_{total}'},'location','northwest') ;
    %grid on ; box on ;
    ylabel('Total Wing Torque [N*m]') ;
    set(gca,'ylim',ylim ) ;
    set(gca,'xlim',xlim ) ;
    %Does it matter that L & R have different flip times???
    plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
    xlabel('Time [ms]')
    title(titleStr) ;
    
s2 = subplot(3,1,2);
    if useOldPlot == 1
        default_ax = gca ;
        set(default_ax,'XTick',[],'YTick',[]) ; %otherwise i get two sets of tick marks
        pos = get(s2,'position') ;
        old_phi_lr = hgload('phi_lr.fig') ;
        ax = gca ;
        set(gca,'fontsize',14) ;
        set(gca,'xlim',xlim) ;
        xlabel('') ;
        title('') ;
        ylabel('Wing Stroke Angle') ;
        ax2 = copyobj(ax, hPhiTorque) ;
        set(ax2,'position',pos) ;
        close(old_phi_lr) ;
    
    elseif useOldPlot == 0
        set(gca,'fontsize',14) ;
        hf = fill(tsfvec , avec,'y') ;
        set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
        hold on
        
        plot(tms,data.phiL,'bo')
        plot(tms,data.phiR,'rx')
        plot(tms,(data.params.beta)*ones(size(tms)),'k-','LineWidth',1.5)
        
        legend({'Pert.','\phi_L','\tau_R','\beta'},'location','southwest') ;
        ylabel('Wing Stroke Angle') ;
        set(gca,'xlim',xlim ) ;
        set(gca,'ylim',[0 200])
        box on; grid on;
        plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
    end

figure(hPhiTorque)
s3 = subplot(3,1,3); %left and right wing torques
    hold on;
    
    htemp1 = plot(tms(i1:i2), pitchTorqueR(i1:i2),'Color',[.5 0 0], 'LineWidth', 2) ;
    htemp2 = plot(tms(i1:i2), pitchTorqueL(i1:i2),'Color',[.5 0 0], 'LineWidth', 2) ;
    ylim2 = get(gca,'ylim') ;
    delete(htemp1) ; delete(htemp2) ; 
    set(gca,'fontsize',14) ;
    
    avec = [ylim2(1) ylim2(1) ylim2(2) ylim2(2) ylim2(1)] ;
    hf = fill(tsfvec , avec,'y') ;
    set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
    
    plot(tms, pitchTorqueR,'Color',[.5 0 0], 'LineWidth', 2) ;
    plot(tms, pitchTorqueL,'Color',[0 0 .5], 'LineWidth', 2) ;
    legend({'Pert.','\tau_R','\tau_L'},'location','northwest') ;
    grid on ; box on ;
    ylabel('Individual Wing Torques [N*m]') ;
    set(gca,'ylim',ylim2 ) ;
    set(gca,'xlim',xlim ) ;
    %Does it matter that L & R have different flip times???
    plotWingstrokeBackground(gca, backFlipTimesR*1000, fwdFlipTimesR*1000, patchColor, true);
%}
%Add stuff to data
%data.pitchTorque = pitchTorque;
%data.pitchTorqueR = pitchTorqueR;
%data.pitchTorqueL = pitchTorqueL;

end