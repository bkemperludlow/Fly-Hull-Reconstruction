%cstart = clock ;
%disp(cstart) ; 

numMovies = 1 ; 

ExprNum = 1 ; %8*ones(1,numMovies) ; %  %[5*ones(1,13) 6*ones(1,13)] ;
MovNum = 5; % [3 5 6 7 8 9 10 11 13 16 18 19 23] ; %%[0 1 2 4 5 7 8 11 12 14 15 17 18 3 4 13 18 20 38 39 40 41 43 44 47 49]  ;
tinArray = nan(1,numMovies) ; %-8*[85 30 50 45 50 35 40 60 80 40 35 55 45 90 85 35 100 70 105 85 70 90 70 85 100 70 ] ; %-8*[85 100 50 50 65 55 70 80 65 40 80 70 50 65] ; %nan(1,17) ;
toutArray = nan(1,numMovies) ; %%8*[75 55 60 80 60 87 58 85 100 90 52 44 72 100 95 39 75 60 95 90 75 60 70 30 58 69] ; %8*[100 100 70 80 60 60 100 60 95 80 80 100 90 100] ; %nan(1,17) ; 
pertFlags = 0 ; %[0 -2 1 2 1 0 1 1 0 2 -2 0 0] ; %%[1 2 1 2 -2 -2 -2 2 -1 1 -2 -2 1 1 -1 2 2 -1 2 2 -1 -1 1 2 2 -2] ;
%-2 = roll left -1 = pitch down, 1 = pitch up , 
%0 = no perturbation, 2 = roll right, 3 = other

TwoFliesXY = zeros(1,numMovies) ;
TwoFliesXZ = zeros(1,numMovies) ;%
TwoFliesYZ = zeros(1,numMovies) ;%

%movieSummaryPath = 'G:\Janelia Flies\kir2.1 flies round 2\' ; 
%[ExprNum, MovNum, tinArray, toutArray, pertFlags, TwoFliesXY, TwoFliesXZ, TwoFliesYZ] = ...
%    generateAnalysisList(movieSummaryPath) ;

%NEED TO SET ROOT PATH PROPERLY!!
rootPath_data = 'H:\Fly Data\Opto Silencing\01_26022018\' ; %'G:\Janelia Flies\kir2.1 flies round 2\' ;
rootPath_save = 'H:\Fly Data\Opto Silencing\01_26022018\' ; %'G:\Janelia Flies\kir2.1 flies round 2\' ;
    
XZ = 2 ;
XY = 3 ;
YZ = 1 ;
runSize = size(ExprNum,2);

for i = 1:runSize %208:248 %1:runSize
    
    if MovNum(i) < 10
        zstr = '00' ;
    elseif MovNum(i) < 100 
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
    cd(rootPath_data)
    ExprNumStr = [zstr2 num2str(ExprNum(i))] ;
    folderPath = dir([ExprNumStr '*']) ;
    
    if length(folderPath) > 1
        disp('Folder name unclear, please check') 
        keyboard ;
    end
   
    dataPath = [rootPath_data folderPath.name] ;
    
    if pertFlags(i) == 1
        savePath = [rootPath_save '\Analysis\Pitch Up'] ; %Don't include Expr/Mov Num
    elseif pertFlags(i) == -1
        savePath = [rootPath_save '\Analysis\Pitch Down'] ;
    elseif pertFlags(i) == 0
        savePath = [rootPath_save '\Analysis\No Perturbation' ] ;
    elseif pertFlags(i) == 2
        savePath = [rootPath_save '\Analysis\Roll Right' ] ;
     elseif pertFlags(i) == -2
        savePath = [rootPath_save '\Analysis\Roll Left' ] ;
    elseif pertFlags(i) == 3
        savePath = [rootPath_save '\Analysis\Unsorted' ] ;
    else
        return ;
    end
    calibrationPath = [dataPath '\calibration'] ;
    
    cinFilenames = cell(1,3) ;
    cinFilenames{XY} = [dataPath strcat('\xy_',movNumStr,'.cin')] ;
    cinFilenames{XZ} = [dataPath strcat('\xz_',movNumStr,'.cin')] ;
    cinFilenames{YZ} = [dataPath strcat('\yz_',movNumStr,'.cin')] ;
    
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
    %if exist(prefixStr,'dir') == 7 % may want to change this back at some point
    %    continue ; 
    %else
    mkdir(prefixStr) ;
    %end
    
    cd(prefixStr) ;
    mkdir('figs');
    DLT_matrix_CSV_filename = [calibrationPath '\calibration_dltCoefs.csv'] ;
    easyWandData_filename = [calibrationPath '\calibration_easyWandData.mat'] ;
    
    format compact
    
    resultsFileName = [prefixStr '_results'] ;
    
    dlt_matrix = load(DLT_matrix_CSV_filename) ; % CSV
    load(easyWandData_filename); % contains easyWandData
    set(0, 'ShowHiddenHandles', 'on')
    close 'easyWand 5'
    
    try
        run beginning
        %run justPitch
    catch
        %keyboard ;
        %try 
        %    run justPitch
        %catch
            cd ..
            diary off
            try
                rmdir(prefixStr,'s')
            catch
                disp(['Could not delete ' prefixStr]) ;
            end
            continue ;
        %end
    end
    %run justPitch
    close all
end

%cend = clock ;
%disp(cend) ; 
