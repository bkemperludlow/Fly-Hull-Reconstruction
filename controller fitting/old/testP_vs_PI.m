rootPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\' ;



logPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Controller Analysis' ;
cd(logPath)

dataLog = importdata('Pitch_Controller_Log.xlsx') ;
ExprNumArray = dataLog.data.Finished(:,1) ;
MovNumArray = dataLog.data.Finished(:,2) ;
pitchTypeArray = dataLog.data.Finished(:,3) ;

ExprNumArray = ExprNumArray(~isnan(ExprNumArray)) ;
MovNumArray = MovNumArray(~isnan(MovNumArray)) ;
pitchTypeArray = pitchTypeArray(~isnan(pitchTypeArray)) ;

Nmovies = length(MovNumArray) ; 
F_test_mat = zeros(Nmovies, 4) ; 

for i = 1:Nmovies
    
    if pitchTypeArray(i) == 1
        dataPath = [rootPath 'Pitch Up'] ;
    elseif pitchTypeArray(i) == -1 ;
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
    elseif MovNum < 100 ;
        zstr = '0' ;
    else
        zstr = '' ;
    end
    
    folderName = ['Expr_' num2str(ExprNum) '_mov_' zstr num2str(MovNum) ] ;
    %fitFilename1 = ['Expr' num2str(ExprNum) 'mov' num2str(MovNum) '_controllerFit_LM'] ;
    fitFilename1 = 'controller_fit_struct_LM.mat' ;
    fitFilename2 = 'controller_fit_struct_LM_Pcontroller.mat' ;
    
    
    cd(folderName)
    
    PI_controller_struct = importdata( fitFilename1 );
    P_controller_struct = importdata( fitFilename2 );
    
    if strcmp(PI_controller_struct.flyType, 'experimental')
        flyType = 1 ; 
    elseif strcmp(PI_controller_struct.flyType, 'control')
        flyType = 2 ; 
    else
        keyboard ; 
    end

    pValue = myFTest(P_controller_struct, PI_controller_struct) ; 
    
    F_test_mat(i,1) = ExprNum ; 
    F_test_mat(i,2) = MovNum ; 
    F_test_mat(i,3) = pValue ; 
    F_test_mat(i,4) = flyType ; 
    
    
end


cd('G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Controller Analysis')
save F_test_mat F_test_mat
