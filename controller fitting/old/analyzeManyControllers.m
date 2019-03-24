%Look at the controller fits for many movies

rootPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\' ;
pertType = 'pitch' ; %or 'roll'

controllerAnalysisStruct = struct('ExprNum', [], 'MovNum', [], 'K_i',[],'K_p',[],'K',[],...
            'deltaT',[],'rms',[],'I_norm',[],'P_norm',[],'medianWingBeat',[],...
            'pertType',[],'flyType',[]);

%for pertType, 1 = pitch up , -1 = pitch down, 2 = roll right, -2 = roll left
% for flyType, 1 = experimental , 2 = control

switch pertType
    case 'pitch'
        logPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Controller Analysis' ;
        cd(logPath)
        
        dataLog = importdata('Pitch_Controller_Log.xlsx') ;
        ExprNumArray = dataLog.data.Finished(:,1) ;
        MovNumArray = dataLog.data.Finished(:,2) ;
        pitchTypeArray = dataLog.data.Finished(:,3) ;
        
        ExprNumArray = ExprNumArray(~isnan(ExprNumArray)) ;
        MovNumArray = MovNumArray(~isnan(MovNumArray)) ;
        pitchTypeArray = pitchTypeArray(~isnan(pitchTypeArray)) ;
        
        
        for i = 1:length(MovNumArray)
            
            if pitchTypeArray(i) == 1
                dataPath = [rootPath 'Pitch Up'] ;
            elseif pitchTypeArray(i) == -1  
                dataPath = [rootPath 'Pitch Down'] ;
            else
                disp('Pitch type not correct. Please check')
                keyboard ; 
            end
            cd(dataPath)
            ExprNum = ExprNumArray(i) ;
            MovNum = MovNumArray(i) ;
            
            if MovNum < 10
                zstr = '00' ;
            elseif MovNum < 100 
                zstr = '0' ;
            else
                zstr = '' ;
            end
            
            folderName = ['Expr_' num2str(ExprNum) '_mov_' zstr num2str(MovNum) ] ;
            %fitFilename1 = ['Expr' num2str(ExprNum) 'mov' num2str(MovNum) '_controllerFit_LM'] ;
            fitFilename1 = 'controller_fit_struct_LM.mat' ;
            fitFilename2 = ['Expr' num2str(ExprNum) 'mov' num2str(MovNum) '_controllerFit'] ;
            dataFilename1 = ['Expr' num2str(ExprNum) 'mov' zstr num2str(MovNum) '_Data_manually_corrected.mat'] ;
            dataFilename2 = ['Expr' num2str(ExprNum) 'mov' zstr num2str(MovNum) '_Data_manually_corrected_onlyPhiFront.mat'] ;
            dataFilename3 = ['Expr_' num2str(ExprNum) '_mov_' zstr num2str(MovNum) '_test.mat'] ;
            
            cd(folderName)
            try 
                load( fitFilename1 )
            catch
                load( fitFilename2 )
            end
            
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
            
            [I_norm, P_norm] = compareControllerContributions(data, controller_fit_struct) ;
            medWingBeat = median(diff(controller_fit_struct.fwdFlipTimes)) ;
            if strcmp({controller_fit_struct.flyType},'control')
                flyType = 2 ;
            elseif strcmp({controller_fit_struct.flyType},'experimental')
                flyType = 1 ;
            else
                flyType = nan ;
            end
            
            controllerAnalysisStruct(i).ExprNum = ExprNum ;
            controllerAnalysisStruct(i).MovNum = MovNum ;
            controllerAnalysisStruct(i).K_i = controller_fit_struct.K_i ;
            controllerAnalysisStruct(i).K_p = controller_fit_struct.K_p ;
            controllerAnalysisStruct(i).K = controller_fit_struct.K ;
            controllerAnalysisStruct(i).deltaT = controller_fit_struct.deltaT ;
            controllerAnalysisStruct(i).rms = controller_fit_struct.rms ;
            controllerAnalysisStruct(i).flyType = flyType ;
            controllerAnalysisStruct(i).I_norm = I_norm ;
            controllerAnalysisStruct(i).P_norm = P_norm ;
            controllerAnalysisStruct(i).medianWingBeat = medWingBeat ;
            controllerAnalysisStruct(i).pertType = pitchTypeArray(i) ;
            controllerAnalysisStruct(i).K_i_CI = controller_fit_struct.K_i_CI ;
            controllerAnalysisStruct(i).K_p_CI = controller_fit_struct.K_p_CI ;
            controllerAnalysisStruct(i).K_CI = controller_fit_struct.K_CI ;
            controllerAnalysisStruct(i).deltaT_CI = controller_fit_struct.deltaT_CI ;
                        
            clear data  
        end
        
        
        cd('G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Controller Analysis')
        save controllerAnalysisStruct_LM controllerAnalysisStruct
        

    case 'roll'
        logPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Roll Controller Analysis' ;
        cd(logPath)
        
        dataLog = importdata('Roll_Controller_Log.xlsx') ;
        ExprNumArray = dataLog.data.Finished(:,1) ;
        MovNumArray = dataLog.data.Finished(:,2) ;
        rollTypeArray = dataLog.data.Finished(:,3) ;
        
        ExprNumArray = ExprNumArray(~isnan(ExprNumArray)) ;
        MovNumArray = MovNumArray(~isnan(MovNumArray)) ;
        rollTypeArray = rollTypeArray(~isnan(rollTypeArray)) ;
        
        for i = 1:length(MovNumArray)
            
            if rollTypeArray(i) == 2
                dataPath = [rootPath 'Roll Right'] ;
            elseif rollTypeArray(i) == -2 ; 
                dataPath = [rootPath 'Roll Left'] ;
            else
                disp('Roll type not correct. Please check')
                keyboard ; 
            end
            cd(dataPath)
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
            dataFilename2 = ['Expr' num2str(ExprNum) 'mov' zstr num2str(MovNum) '_Data_manually_corrected_onlyPhi.mat'] ;
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
            
            [I_norm, P_norm] = compareControllerContributions(data, controller_fit_struct) ;
            medWingBeat = 2*median(diff(controller_fit_struct.phiAmpTimes)) ;
            if strcmp({controller_fit_struct.flyType},'control')
                flyType = 2 ;
            elseif strcmp({controller_fit_struct.flyType},'experimental')
                flyType = 1 ;
            else
                flyType = nan ;
            end
            
            controllerAnalysisStruct(i).ExprNum = ExprNum ;
            controllerAnalysisStruct(i).MovNum = MovNum ;
            controllerAnalysisStruct(i).K_i = controller_fit_struct.K_i ;
            controllerAnalysisStruct(i).K_p = controller_fit_struct.K_p ;
            %controllerAnalysisStruct(i).K = controller_fit_struct.K ;
            controllerAnalysisStruct(i).deltaT = controller_fit_struct.deltaT ;
            controllerAnalysisStruct(i).rms = controller_fit_struct.rms ;
            controllerAnalysisStruct(i).flyType = flyType ;
            controllerAnalysisStruct(i).I_norm = I_norm ;
            controllerAnalysisStruct(i).P_norm = P_norm ;
            controllerAnalysisStruct(i).medianWingBeat = medWingBeat ;
            controllerAnalysisStruct(i).pertType = rollTypeArray(i) ;
                        
            clear data  
        end
         %{
        cd('G:\Janelia Flies\kir2.1 flies round 2\Analysis\Roll Controller Analysis')
        save controllerAnalysisStruct controllerAnalysisStruct
        %}
        
end


