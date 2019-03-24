%mk2: Trying to figure out how to automate categorization 

%cstart = clock ;
%disp(cstart) ; 

%testing running analyses automatically via file detection
%changed to function so watchForVids can call it
%mk2: Changed so pertFlags no longer needed; analysis will guess
%perturbation type and place it into corresponding folder
function runManyAnalysesAuto_forClust(MovNum, ExprNum, pathStruct)

%ExprNum =  14 ; %[5*ones(1,13) 6*ones(1,13)] ;
%MovNum =  1 ;%[0 1 2 4 5 7 8 11 12 14 15 17 18 3 4 13 18 20 38 39 40 41 43 44 47 49]  ;
tinArray = nan ; %-8*[85 30 50 45 50 35 40 60 80 40 35 55 45 90 85 35 100 70 105 85 70 90 70 85 100 70 ] ; %-8*[85 100 50 50 65 55 70 80 65 40 80 70 50 65] ; %nan(1,17) ;
toutArray = nan ; %8*[75 55 60 80 60 87 58 85 100 90 52 44 72 100 95 39 75 60 95 90 75 60 70 30 58 69] ; %8*[100 100 70 80 60 60 100 60 95 80 80 100 90 100] ; %nan(1,17) ; 
%pertFlags = 0 ; %[1 2 1 2 -2 -2 -2 2 -1 1 -2 -2 1 1 -1 2 2 -1 2 2 -1 -1 1 2 2 -2] ;
%-2 = roll left -1 = pitch down, 1 = pitch up , 
%0 = no perturbation, 2 = roll right, 3 = other
rootPath = pathStruct.root ; 

falseTriggerFlag = false ;
errorFlag = false ;

TwoFliesXY = 0 ;
TwoFliesXZ = 0 ;%
TwoFliesYZ = 0 ;%

%NEED TO SET ROOT PATH PROPERLY!!
%rootPath = 'D:\Raymond Analysis Test\Raw Movies\' ;
%rootPath = 'D:\Raymond Analysis Test\test automation\' ;
%rootPath = 'D:\Raymond Analysis Test\camera test automation\' ;

XZ = 2 ;
XY = 3 ;
YZ = 1 ;
runSize = size(ExprNum,2);

for i = 1:runSize
    
    if MovNum(i) < 10
        zstr = '00' ;
    elseif MovNum(i) < 100 ;
        zstr = '0' ;
    else
        zstr = '' ;
    end
    movNumStr = [zstr num2str(MovNum(i))] ;
    
    if ExprNum(i) < 10
        zstr2 = '0' ;
    else
        zstr2 = '' ;
    end
    cd(rootPath)
    ExprNumStr = [zstr2 num2str(ExprNum(i))] ;
    folderPath = dir([ExprNumStr '*']) ;
    
    if length(folderPath) > 1
        disp('Folder name unclear, please check') 
        keyboard ;
    end
   
    %dataPath = [rootPath folderPath.name] ;
    dataPath = rootPath ; 
    
%     if pertFlags(i) == 1
%         savePath = [rootPath '\Analysis\Pitch Up'] ; %Don't include Expr/Mov Num
%     elseif pertFlags(i) == -1
%         savePath = [rootPath '\Analysis\Pitch Down'] ;
%     elseif pertFlags(i) == 0
%         savePath = [rootPath '\Analysis\No Perturbation' ] ;
%     elseif pertFlags(i) == 2
%         savePath = [rootPath '\Analysis\Roll Right' ] ;
%      elseif pertFlags(i) == -2
%         savePath = [rootPath '\Analysis\Roll Left' ] ;
%     elseif pertFlags(i) == 3
%         savePath = [rootPath '\Analysis\Other' ] ;
%     else
%         return ;
%     end
    savePath = pathStruct.save ;
    calibrationPath = pathStruct.calibration ;
    errorPath = [rootPath 'errorlog.txt'] ;
    possibleFTPath = pathStruct.possibleFT ;
    
        
    cinFilenames = cell(1,3) ;
    cinFilenames{XY} = [dataPath strcat('\xy_',movNumStr,'.cin')] ;
    cinFilenames{XZ} = [dataPath strcat('\xz_',movNumStr,'.cin')] ;
    cinFilenames{YZ} = [dataPath strcat('\yz_',movNumStr,'.cin')] ;
    xmlFilenames = cell(1,3) ;
    xmlFilenames{XY} = [dataPath strcat('\xy_',movNumStr,'.xml')] ;
    xmlFilenames{XZ} = [dataPath strcat('\xz_',movNumStr,'.xml')] ;
    xmlFilenames{YZ} = [dataPath strcat('\yz_',movNumStr,'.xml')] ;
    
    expr_number = ExprNum(i) ;
    twoFlies_xy = TwoFliesXY(i) ;
    twoFlies_xz = TwoFliesXZ(i) ;
    twoFlies_yz = TwoFliesYZ(i) ;
    
    % find common tin and tout
    tin = tinArray(i) ; % max(allTin) ;
    tout = toutArray(i); %-611 ;% 899 % min(allTout) ;
    
    %dataPath = 'F:\Pitch\06_220514\' ;
    %dataPath = 'F:\Pitch\07_230514\' ;
    
    prefixStr = ['Expr_' num2str(ExprNum(i)) '_mov_' movNumStr ];
    cd(savePath);
    mkdir(prefixStr) ;
    
    cd(prefixStr) ;
    mkdir('figs');
    DLT_matrix_CSV_filename = [calibrationPath '\calibration_dltCoefs.csv'] ;
    easyWandData_filename = [calibrationPath '\calibration_easyWandData.mat'] ;
    
    %format compact
    
    resultsFileName = [prefixStr '_results'] ;
    
    dlt_matrix = load(DLT_matrix_CSV_filename) ; % CSV
    load(easyWandData_filename); % contains easyWandData
    set(0, 'ShowHiddenHandles', 'on')
    close 'easyWand 5'
    
    try
        %run beginning_mk2
        %run beginning_mk3
        run beginning_forClust
        if falseTriggerFlag
            %moves false trigger files into corresponding analysis folder
            movefile(strcat(savePath, '\', prefixStr), possibleFTPath) ;
            movefile(cinFilenames{XY}, [possibleFTPath '\' prefixStr]) ;
            movefile(cinFilenames{XZ}, [possibleFTPath '\' prefixStr]) ;
            movefile(cinFilenames{YZ}, [possibleFTPath '\' prefixStr]) ;
            movefile(xmlFilenames{XY}, [possibleFTPath '\' prefixStr]) ;
            movefile(xmlFilenames{XZ}, [possibleFTPath '\' prefixStr]) ;
            movefile(xmlFilenames{YZ}, [possibleFTPath '\' prefixStr]) ;
            close all ;
            return ; %continue
        elseif errorFlag
            return ; %continue
        else
            data = batchCalcAngles_freeFlight_mk2(ExprNum, movNumStr, savePath) ;
        end
        %run justPitch
    catch exception
        %keyboard ;
        %try 
        %    run justPitch
        %catch
            disp(getReport(exception, 'basic'))
            cd ..
            diary off
            try
                rmdir(prefixStr,'s')
            catch
                disp(['Could not delete ' prefixStr]) ;
            end
            return ; %continue
        %end
    end
    %run justPitch
    close all
    %tries to guess perturbation type and moves it to corresponding folder
    try 
        run guessPerturbation_mk2
    catch 
        pertStr = 'Error' ; 
    end
    pertLogChange = strcat(num2str(MovNum), pertStr) ;  
    fileID = fopen(strcat(savePath,'\','pertLog.txt'),'a+') ;
    fprintf(fileID, '%s\r\n', pertLogChange) ;
    fclose(fileID) ;
    movefile(strcat(savePath,'\',prefixStr), strcat(rootPath,'Analysis\Unsorted\')) ;
end
%cend = clock ;
%disp(cend) ; 
