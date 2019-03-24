rootPath = 'D:\Box Sync\VNC Motor Lines\23_21062018\Analysis\Unsorted' ;
highestNumFolder = 41;
ExprNum = 23 ;
pertType = nan ;
flyType = 2 ;
plotFlag1 = false ;
plotFlag2 = true ;
debugFlag1 = false ;
debugFlag2 = false ;
saveFlag = false ;
smoothFlag = true ;

defineConstantsScript ;

for MovNum = 1:highestNumFolder
    try   
        %load data structure
        data = loadPertDataStruct(rootPath, ExprNum,MovNum,pertType);
        [deltaAlphaMean, c_yaw, midWingBeatTimes, deltaAlphaSEM] = ...
            getYawControllerData(data, debugFlag1) ;
        
        %get body yaw raw data
        bodyYaw = data.anglesLabFrame(:,PHIB) ;
        
        %get time 
        fps = data.params.fps ;
        startTime = data.params.startTrackingTime ;
        endTime   = data.params.endTrackingTime ;
        Np = data.Nimages ;
        t = (startTime:endTime) / fps  ; % in SEC
        
        %get correct time range for the model yaw
        t1_idx = find(t < midWingBeatTimes(1),1,'last') ;
        t2_idx = find(t > midWingBeatTimes(end),1,'first') ;
        %t0_idx = find(t == 0) ;
        
        if smoothFlag
            deltaAlphaMean = smooth(deltaAlphaMean,3) ;
        end
        
        %plot attack angle difference and model yaw
        yawplot = figure;
        subplot(2,1,1);
        
        %left y axis
        yyaxis left;
        errorbar(1000*midWingBeatTimes, deltaAlphaMean,deltaAlphaSEM,'ko-');
        ylabel('\Delta\alpha');
        axis tight;
       
        %right y axis
        yyaxis right;
        plot(1000*t(t1_idx:t2_idx), c_yaw(t(t1_idx:t2_idx)) - c_yaw(0));
        ylabel('Model Yaw');
        title(['Model Yaw and Attack Angle Difference ' num2str(MovNum, '%03d') ]);
        axis tight;
        xlabel('Time');
        
%{
        yyaxis right; ylimr = get(gca,'Ylim');ratio = ylimr(1)/ylimr(2);
        yyaxis left; yliml = get(gca,'Ylim');
        if yliml(2)*ratio<yliml(1)
            set(gca,'Ylim',[yliml(2)*ratio yliml(2)])
        else
            set(gca,'Ylim',[yliml(1) yliml(1)/ratio])
        
        end
        %}
 
       
        %plot raw data
        subplot(2,1,2);
        plot(1000*t(t1_idx:t2_idx), bodyYaw(t1_idx:t2_idx), '.');
        title('Raw Data');
        ylabel('Yaw');
        xlabel('Time');
        axis tight;
        folderName = ['Expr_' num2str(ExprNum) '_mov_' num2str(MovNum, '%03d') ] ;
        dataPath = [rootPath '\' folderName];
        saveas(yawplot,[dataPath '\' 'Control_Model'],'png') ; 
        saveas(yawplot,[rootPath '\' 'Control_Model' num2str(MovNum, '%03d')]) 
    catch
        disp(['wrong' num2str(MovNum, '%03d')]);
    end
end
%{
for folder = dir'
    display(folder.name);
    dataPath = [rootPath '\' folder.name];
    display(dataPath)
    dataFileName1 = [dataPath '\' folder.name '_Data_manually_corrected.mat'];
    dataFileName2 = [dataPath '\' folder.name '_test.mat'];
    display(dataFileName1);
    display(dataFileName2);
    if exist(dataFileName1)
        importdata(dataFileName1);
        figure;
        yyaxis left;
        errorbar(1000*midWingBeatTimes, deltaAlphaMean,deltaAlphaSEM,'ko-');
        yyaxis right;
        plot(1000*t(t_start_ind:end), c_yaw(t(t_start_ind:end)));
    elseif exist(dataFileName2)
        load(dataFileName2);
        figure;
        yyaxis left;
        errorbar(1000*midWingBeatTimes, deltaAlphaMean,deltaAlphaSEM,'ko-');
        yyaxis right;
        plot(1000*t(t_start_ind:end), c_yaw(t(t_start_ind:end)));
    end
    display('end')
    end
%}

