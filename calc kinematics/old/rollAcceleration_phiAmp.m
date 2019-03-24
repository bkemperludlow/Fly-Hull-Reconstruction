rootPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\' ;
%cd(rootPath)

logPath = [rootPath 'Roll Controller Analysis'] ;
cd(logPath)
dataLog = importdata('Roll_Controller_Log.xlsx') ;

ExprNumArray = dataLog.data.Finished(:,1) ;
MovNumArray = dataLog.data.Finished(:,2) ;
pertTypeArray = dataLog.data.Finished(:,3) ;

ExprNumArray = ExprNumArray(~isnan(ExprNumArray)) ; 
MovNumArray = MovNumArray(~isnan(MovNumArray)) ; 
pertTypeArray = pertTypeArray(~isnan(pertTypeArray)) ; 

defineConstantsScript

accelerationAndPhiAmpStruct = struct('ExprNum', [], 'MovNum', [], 't',[],'bodyRoll',[],...
    'rollVelocity',[],'rollAcceleration',[],'phiAmpTimes', [], 'phiAmpDiff',[],...
    'maxPhiAmpDiff', [],'maxRollAccel',[],'deltaRoll',[],'peakInd',[],'pertType',[],'flyType',[]);

%d1 = designfilt('lowpassiir','FilterOrder',8,'SampleRate',8000, ...
%        'HalfPowerFrequency',100,'DesignMethod','butter'); %hpf = 100
    
%cc = 1 ; 

for i = 1:length(MovNumArray)
    % load data 
    if pertTypeArray(i) == 2
        cd([rootPath 'Roll Right']) 
    elseif pertTypeArray(i) == -2
        cd([rootPath 'Roll Left']) 
    else 
        keyboard ; 
    end 
    
    ExprNum = ExprNumArray(i) ;
    MovNum = MovNumArray(i) ;

    if MovNum < 10
        zstr = '00' ;
    elseif MovNum < 100 ;
        zstr = '0' ;
    else
        zstr = '' ;
    end

    folderName = ['Expr_' num2str(ExprNum) '_mov_' zstr num2str(MovNum) ] ;
    fitFilename = ['Expr' num2str(ExprNum) 'mov' num2str(MovNum) '_controllerFit'] ;
    dataFilename1 = ['Expr' num2str(ExprNum) 'mov' zstr num2str(MovNum) '_Data_manually_corrected.mat'] ;
    dataFilename2 = ['Expr' num2str(ExprNum) 'mov' zstr num2str(MovNum) '_Data_manually_corrected_onlyPhiFront.mat'] ;
    dataFilename3 = ['Expr_' num2str(ExprNum) '_mov_' zstr num2str(MovNum) '_test.mat'] ;
    
    cd(folderName) 
    load( fitFilename )
    if exist(dataFilename1,'file') == 2
        load(dataFilename1)
    elseif exist(dataFilename2,'file') == 2
        load(dataFilename2)
    elseif exist(dataFilename3,'file') == 2
        load(dataFilename3)
    else
        disp('No data file')
        continue ;
    end
    
    if strcmp({controller_fit_struct.flyType},'control')
        flyType = 2 ;
    elseif strcmp({controller_fit_struct.flyType},'experimental') 
        flyType = 1 ;
    else 
        flyType = nan ;
    end
    
    %t = controller_fit_struct.t ;
    t = (data.params.startTrackingTime : data.params.endTrackingTime ) / 8000 ; 
    phiAmpTimes = controller_fit_struct.phiAmpTimes ;
    phiAmpDiff = controller_fit_struct.phiAmpDiff ;
    
    %get body roll info
    
    sp_rho = controller_fit_struct.sp_rho ; 
    bodyRoll = fnval(sp_rho, t) ;
    rollVel = fnval(fnder(sp_rho,1), t) ;
    rollAccel = fnval(fnder(sp_rho,2), t) ;
    bodyRoll_prePert = mean(bodyRoll(t <= 0 & t > -0.01)) ; %pitch_temp(t0_ind) ;
    deltaRoll = abs(bodyRoll - bodyRoll_prePert) ; 
    [pks, locs] = findpeaks(deltaRoll(t > 0.007)) ;
    postPulse_ind = find(t  > 0.007, 1, 'first') ;
    
    %if numel(pks) > 1
    %    if (pks(2)/pks(1) > 3) && (locs(2) < 180)
    %        whichPeak = 2 ;
    %    else
    %        whichPeak = 1 ;
    %    end
    %else
    %    whichPeak = 1 ; 
    %end
    ind_peak = locs(1) + postPulse_ind ;
    
    %check that the sign is right
    deltaRoll_signed = bodyRoll(ind_peak) - bodyRoll_prePert ; 
    if sign(deltaRoll_signed) ~= sign( pertTypeArray(i) ) 
        disp('Peak does not match pert type')
        disp(ExprNum)
        disp(MovNum)
        %keyboard ;
    end
    t_peak = t(ind_peak) ; 
    maxRollAccel = rollAccel(ind_peak) ; 
    
    %get changes in wing kinematics
    [~,maxPhiAmpDiffInd] = max(abs(phiAmpDiff)) ; 
    maxPhiAmpDiff = phiAmpDiff(maxPhiAmpDiffInd) ; 
    
    if (0)
        t_min = 1000*phiAmpTimes(1) ;
        t_max = 1000*phiAmpTimes(end) ;
        
        figure ;
        subplot(2,2,1)
        hold on
        plot(1000*t, bodyRoll,'b-','LineWidth',1.5)
        %plot(1000*t, bodyPitch, 'k.')
        plot(1000*t_peak, bodyRoll(ind_peak), 'ro', 'MarkerFaceColor','r')
        ylabel('\rho [deg]')
        
        subplot(2,2,3)
        plot(1000*t, rollVel,'m-','LineWidth',1.5)
        ylabel('Roll Vel [deg/s]')
        
        subplot(2,2,2)
        plot(1000*t, rollAccel,'r-','LineWidth',1.5)
        hold on
        plot(1000*t_peak, rollAccel(ind_peak),'bx')
        set(gca,'xlim',[t_min t_max])
        ylabel('Roll Accel [deg/s^2]')
        xlabel('Time [ms]')
        
        subplot(2,2,4)
        hold on
        plot(1000*phiAmpTimes, phiAmpDiff, 'ko-')
        %plot(1000*fwdFlipTimes, phiFront,'kv','MarkerFaceColor', 'k')
        %plot(1000*backFlipTimes, phiBack,'k^','MarkerFaceColor', 'k')
        set(gca,'xlim',[t_min t_max])
        ylabel('\Phi_{diff} [deg]')
        
    end
    
    accelerationAndPhiAmpStruct(i).ExprNum = ExprNum ;
    accelerationAndPhiAmpStruct(i).MovNum = MovNum ;
    accelerationAndPhiAmpStruct(i).t = t ;
    accelerationAndPhiAmpStruct(i).bodyRoll = bodyRoll ;
    accelerationAndPhiAmpStruct(i).deltaRoll = deltaRoll ;
    accelerationAndPhiAmpStruct(i).rollVelocity = rollVel ;
    accelerationAndPhiAmpStruct(i).rollAcceleration = rollAccel ;
    accelerationAndPhiAmpStruct(i).phiAmpTimes = phiAmpTimes ;
    accelerationAndPhiAmpStruct(i).phiAmpDiff = phiAmpDiff ;
    accelerationAndPhiAmpStruct(i).maxPhiAmpDiff = maxPhiAmpDiff ;
    accelerationAndPhiAmpStruct(i).maxRollAccel = maxRollAccel ;
    accelerationAndPhiAmpStruct(i).peakInd = ind_peak ;
    accelerationAndPhiAmpStruct(i).flyType = flyType ;
    accelerationAndPhiAmpStruct(i).pertType = pertTypeArray(i) ; 
    
    %cc = cc + 1 ;
    clear data 
end



%{
cd(logPath)
save accelerationAndPhiAmpStruct accelerationAndPhiAmpStruct
%}