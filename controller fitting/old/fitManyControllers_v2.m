
logPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\Pitch Controller Analysis' ;
cd(logPath)
rootPath = 'G:\Janelia Flies\kir2.1 flies round 2\Analysis\' ;

dataLog = importdata('Pitch_Controller_Log.xlsx') ;
ExprNumArray = dataLog.data.Finished(:,1) ;
MovNumArray = dataLog.data.Finished(:,2) ;
flyTypeArray = dataLog.data.Finished(:,5) ;
pitchTypeArray = dataLog.data.Finished(:,3) ;

ExprNumArray = ExprNumArray(~isnan(ExprNumArray)) ;
MovNumArray = MovNumArray(~isnan(MovNumArray)) ;
flyTypeArray = flyTypeArray(~isnan(flyTypeArray)) ;
pitchTypeArray = pitchTypeArray(~isnan(pitchTypeArray)) ;

for i = 1:length(MovNumArray)
    ExprNum = ExprNumArray(i) ;
    MovNum = MovNumArray(i) ;
    flyType = flyTypeArray(i) ;
    
    if pitchTypeArray(i) == 1
        dataPath = [rootPath 'Pitch Up'] ;
    elseif pitchTypeArray(i) == -1 
        dataPath = [rootPath 'Pitch Down'] ;
    else
        disp('Pitch type not correct. Please check')
        keyboard ;
    end
    cd(dataPath)
    %ExprNum = ExprNumArray(i) ;
    %MovNum = MovNumArray(i) ;
    
    if MovNum < 10
        zstr = '00' ;
    elseif MovNum < 100 
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
    %ExprNum = 2 ;
    %MovNum = 120;
    %flyType = 2 ; %1 = experimental , 2 = control ;
    debugFlag = true ;
    plotFlag = true ;
    numInitializations = 100 ;
    
    lowerBound = [-3 0 .004 -4 ] ; %[K_i,  K_p,  deltaT,  K]
    upperBound =  [3 .3 .025 4 ] ;
    
    if numInitializations <= 1
        %as intial guesses, use more or less what's in the pitch paper
        K_i_guess = 0.0 ;%0.3 ; %0.3 ; %unitless
        K_p_guess = 0.008 ; %seconds
        deltaT_guess = 0.006 ; %seconds
        K_guess = 0 ; %degrees
        
        %x_0 = [K_i_guess K_p_guess deltaT_guess K_guess deltaT_2_guess]; %initial guess for x
        initialGuess = [K_i_guess K_p_guess deltaT_guess K_guess ]; %initial guess for x
        
        controller_fit_struct = controllerFit_v3(data, ExprNum, MovNum, ...
            flyType, initialGuess, lowerBound, upperBound, debugFlag, plotFlag) ;
        
    else
        fit_struct_cell = cell(numInitializations,1) ;
        rmseArray = nan(numInitializations,1) ;
        K_i_Array = nan(numInitializations,1) ;
        for j = 1:numInitializations
            K_i_guess = (upperBound(1) - lowerBound(1))*rand(1) + lowerBound(1) ;%0.3 ; %0.3 ; %unitless
            K_p_guess = (upperBound(2) - lowerBound(2))*rand(1) + lowerBound(2) ; %seconds
            deltaT_guess = (upperBound(3) - lowerBound(3))*rand(1) + lowerBound(3) ; %seconds
            K_guess = (upperBound(4) - lowerBound(4))*rand(1) + lowerBound(4) ; %degrees
            
            %x_0 = [K_i_guess K_p_guess deltaT_guess K_guess deltaT_2_guess]; %initial guess for x
            initialGuess = [K_i_guess K_p_guess deltaT_guess K_guess ]; %initial guess for x
            
            controller_fit_struct = controllerFit_v3(data, ExprNum, MovNum, ...
                flyType, initialGuess, lowerBound, upperBound, false, false) ;
            
            fit_struct_cell{j} = controller_fit_struct ;
            rmseArray(j) = controller_fit_struct.rms ;
            K_i_Array(j) = controller_fit_struct.K_i ;
        end
        
        [~, bestFitInd] = min(rmseArray) ;
        controller_fit_struct = fit_struct_cell{bestFitInd} ;
        
        if plotFlag
            t = controller_fit_struct.t/1000 ;
            manualCorrRange = t ;
            c_pitch = controller_fit_struct.c_pitch ;
            deltaT = controller_fit_struct.deltaT ;
            K_i = controller_fit_struct.K_i ;
            K_p = controller_fit_struct.K_p ;
            K = controller_fit_struct.K ;
            phiFront_meansub = controller_fit_struct.deltaPhiFront ;
            fwdFlipTimes = controller_fit_struct.fwdFlipTimes ;
            backFlipTimes = controller_fit_struct.backFlipTimes ;
            
            pitchVel = differentiate(c_pitch,t-deltaT) ;
            controlPred = K + K_p*pitchVel...
                + K_i*(c_pitch(t-deltaT)- c_pitch(0)) ;
            
            proportionalTerm = K_p*pitchVel ;
            integralTerm = K_i*(c_pitch(t-deltaT)- c_pitch(0)) ;
            
            ylim = [min(phiFront_meansub)-3 , max(phiFront_meansub)+3] ;
            tsfvec = [0 7 7 0 0] ;
            patchColor = [1 1 1 ] * 0.8;
            
            if flyType == 1
                plotColor = [.7 0 0 ] ;
            elseif flyType == 2
                plotColor = [0 .7 0 ] ;
            else
                plotColor = .5*[1 1 1] ;
            end
            
            hcontrib = figure ;
            hold on
            %set(gcf, 'PaperPositionMode', 'manual');
            set(gcf, 'Position', [500 500 420 140]);
            set(gcf,'PaperPositionMode','auto')
            %set(hcontrib,'name',['Expr ' num2str(ExprNum) ' Mov ' num2str(MovNum) ' Contributions'],'numbertitle','off')
            
            avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1)] ;
            hf = fill(tsfvec , avec,'y') ;
            set(hf,'facecolor',[255 238 170]/255,'facealpha',1,'edgecolor','none') ;
            set(hf,'HandleVisibility','off')
            
            plot(1000*t,controlPred,'Color',plotColor,'LineWidth',2.5)
            plot(1000*t,proportionalTerm,'Color',0.6*[1 1 1],'LineWidth',1.5)
            plot(1000*t,integralTerm,'k--','LineWidth',1.5)
            %shadedErrorBar(t*1000, controlPred, 2*ones(size(controlPred)),{'-','LineWidth',2.5,'Color',plotColor},1) ;
            %errorbar(1000*fwdFlipTimes, phiFront_meansub, 2*ones(size(phiFront_meansub)), 'ko','markerfacecolor','k') ;
            errorbar(1000*fwdFlipTimes, phiFront_meansub, 2*ones(size(phiFront_meansub)),'ko','markerfacecolor','k') ;
            plotWingstrokeBackground(gca, backFlipTimes*1000, fwdFlipTimes*1000, patchColor, true);
            axis tight ;
            set(gca, 'xlim', 1000*[manualCorrRange(1), max(fwdFlipTimes)])
            set(gca,'ylim',ylim)
            xlabel('Time [ms]')
            ylabel('\Delta \Phi_{front} [deg]')
            %legend({'PI','P','I','Data'})
            
            figure ;
            plot(K_i_Array,rmseArray, 'ko','MarkerFaceColor',plotColor)
            
        end
        
    end
    
    disp('K_i:') ; disp(controller_fit_struct.K_i) ;
    disp('K_p:') ; disp(controller_fit_struct.K_p) ;
    disp('deltaT:') ; disp(controller_fit_struct.deltaT) ;
    disp('K:') ; disp(controller_fit_struct.K) ;
    disp('RMSE:') ; disp(controller_fit_struct.rms) ;
    
    filename = ['Expr' num2str(ExprNum) 'mov' num2str(MovNum) '_controllerFit_NEW'] ;
    
    
    %save(filename, 'controller_fit_struct') ;
    %savefig(gcf,[filename '.fig']) ;
    
end

%--------------------------------------------
%{
rootPath = 'G:\Janelia Flies\kir2.1 flies\Analysis\Pitch Up' ;
cd(rootPath)

%pitchUpControllerFits = struct('ExprNum',[],'MovNum',[],'K_i',[],'K_p',[],...
%    'deltaT',[],'K',[],'rms',[],'sp_pitch',[],'pitchEstErr',[],'sp_phiR',[],...
%    'sp_phiL',[],'phiEstErr',[],'deltaPhiFront',[],'fwdFlipTimes',[],'t',[],...
%    'flyType',[]) ;

ExprNumArray = [4 4 4 4 4 11 12 12 3 3 3 10 10 13] ;
MovNumArray = [4 10 137 147 202 20 11 64 10 30 37 45 38 134 ] ;

%controlInd = find(strcmp({pitchUpControllerFits(:).flyType},'control')) ;
K_i_control = [] ;
K_p_control = [] ;
deltaT_control = [] ;
rms_control = [] ;

%experimentalInd = find(strcmp({pitchUpControllerFits(:).flyType},'experimental')) ;
%experimentalInd = 6:8 ;
K_i_experimental = [] ;
K_p_experimental = [] ;
deltaT_experimental = [] ;
rms_experimental = [] ;


for j = 1:length(MovNumArray)
    cd(rootPath)
    ExprNum = ExprNumArray(j) ;
    MovNum = MovNumArray(j) ;

    if MovNum < 10
        zstr = '00' ;
    elseif MovNum < 100 ;
        zstr = '0' ;
    else
        zstr = '' ;
    end

    folderName = ['Expr_' num2str(ExprNum) '_mov_' zstr num2str(MovNum) ] ;
    filename = ['Expr' num2str(ExprNum) 'mov' num2str(MovNum) '_controllerFit'] ;
    
    cd(folderName)
    load( filename )
    
    if strcmp({controller_fit_struct.flyType},'control')
        K_i_control = [ K_i_control, controller_fit_struct.K_i ] ;
        K_p_control = [ K_p_control, controller_fit_struct.K_p ] ;
        deltaT_control = [ deltaT_control, controller_fit_struct.deltaT ] ;
        rms_control = [ rms_control, controller_fit_struct.rms ] ;
    elseif strcmp({controller_fit_struct.flyType},'experimental')
        K_i_experimental = [ K_i_experimental, max([controller_fit_struct.K_i, 0]) ] ;
        K_p_experimental = [ K_p_experimental, controller_fit_struct.K_p ] ;
        deltaT_experimental = [ deltaT_experimental, controller_fit_struct.deltaT ] ;
        rms_experimental = [ rms_experimental, controller_fit_struct.rms ] ;
    end

end


figure ; hold on

scatter3(K_i_control, 1000*K_p_control, 1000*deltaT_control, 20*rms_control , ...
    [0 0.7 0], 'filled')
scatter3(K_i_experimental, 1000*K_p_experimental, 1000*deltaT_experimental, 20*rms_experimental , ...
    [0.7 0 0], 'filled')

axis tight
grid on
box on

xlabel('K_i [unitless]')
ylabel('K_p [ms]')
zlabel('\Delta T [ms]')
set(gca, 'fontsize',14)

figure ; hold on

scatter(K_i_control, 1000*K_p_control, 20*rms_control , ...
    [0 0.7 0], 'filled')
scatter(K_i_experimental, 1000*K_p_experimental, 20*rms_experimental , ...
    [0.7 0 0], 'filled')

axis tight
grid on
%box on

xlabel('K_i [unitless]')
ylabel('K_p [ms]')
%zlabel('\Delta T [ms]')
set(gca, 'fontsize',14)
%}

