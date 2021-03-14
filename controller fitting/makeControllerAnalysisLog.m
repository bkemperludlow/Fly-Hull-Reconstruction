% -------------------------------------------------------------------------
% Script to automatically generate a log file (xlsx) for all controller
% fits of movies from experiments in a specified directory
% -------------------------------------------------------------------------
% ------------------------
%% path info and params
rootPath = 'D:\Fly Data\VNC Motor Lines\'  ; % 'D:\Fly Data\VNC Motor Lines\' % 'D:\Fly Data\Janelia Flies\kir round 2\'
pertTypes = {'pitchDown', 'pitchUp', 'rollLeft', 'rollRight', 'noPert'} ;
pertTypeNames = {'Pitch Down', 'Pitch Up', 'Roll Left', 'Roll Right', ...
    'No Pert'} ;
saveFilename = fullfile(rootPath, 'AnalysisLogAuto.xlsx') ;
overWriteFlag = true ;

if exist(saveFilename,'file') && ~overWriteFlag
    disp('Already created this log file')
    return
end
% ------------------------------------------------------
%% get directory of all experiment folders in rootPath
exprDir = dir(rootPath) ;
exprDir = exprDir(3:end) ; % get rid of '..' and '.'
exprDir = exprDir([exprDir(:).isdir]) ; % take only folders

% take only directories of the form \d\d_*
exprNumbers = arrayfun(@(x) str2double(x.name(1:2)), exprDir) ;
exprIdx = ~isnan(exprNumbers) ;
exprDir = exprDir(exprIdx) ;
exprNumbers = exprNumbers(exprIdx) ;

% ------------------------------------------------------
%% load catalog file that contains driver info
catalogDir = dir(fullfile(rootPath, '*catalog.csv')) ; 
if (length(catalogDir) == 1)
    catalogFlag = true ; 
    % load in catalog data
    catalog_fn = fullfile(catalogDir(1).folder, catalogDir(1).name) ; 
    opts = detectImportOptions(catalog_fn);
    catalogData = readtable(catalog_fn, opts) ;
    catalogData = table2struct(catalogData) ; 
else
    fprintf('Could not find catalog file for %s -- skipping \n', rootPath)
    catalogFlag = false ; 
end
% -------------------------------------------------------------------------
%% loop through experiment folders and tabulate the controller fit results
warning( 'off', 'MATLAB:xlswrite:AddSheet' ) ; % turn off warning
for i = 1:length(exprDir)
    % get experiment number and path
    exprPath = fullfile(exprDir(i).folder, exprDir(i).name) ;
    pathStruct = generatePathStruct(exprPath) ;
    ExprNum = pathStruct.ExprNum ;
    
    % initialize structure to store analysis results for this experiment
    analysisLogStruct = struct() ;
    cc = 1 ; % counter for filling in struct
    
    % ----------------------------------------------------------
    % use catalog genotype information, if available
    if catalogFlag
        c_ind = find([catalogData.ExprNum] == ExprNum) ;
        if (length(c_ind) ~= 1)
            driver = [] ;
            effector = [] ;
            mn_name = [] ;
        else
            driver = catalogData(c_ind).Driver ;
            effector = catalogData(c_ind).Effector ;
            mn_name = catalogData(c_ind).MN(2:(end-1)) ;
        end
    end
    
    % -------------------------------------------------------------------
    % loop through pert types and check movie folders that correspond to
    % this pert type (i'm assuming that these are sorted by pert type)
    for j = 1:length(pertTypes)
       pt = pertTypes{j} ; 
       pertPath = pathStruct.(pt) ; 
       pt_name = pertTypeNames{j} ; 
       
       % get directory for perturbation folder (all movies) and controller
       % fits. Then we'll record all of these
       pertDir = dir(fullfile(pertPath, 'Expr*')) ; 
       controlFitDir = dir(fullfile(pertPath, '**', 'controller_fit_struct*')) ; 
       
       allMovNums = arrayfun(@(x) str2double(x.name(end-2:end)),pertDir) ;
       fitMovNums = arrayfun(@(x) str2double(x.folder(end-2:end)),controlFitDir) ;
       
       [~, match_ind, ~] = intersect(allMovNums, fitMovNums) ; 
       
       % loop through each movie
       for k = 1:length(allMovNums)
          analysisLogStruct(cc).ExprNum = ExprNum ; 
          analysisLogStruct(cc).MovNum = allMovNums(k) ; 
          analysisLogStruct(cc).pertType = pt_name ; 
          analysisLogStruct(cc).controlFitBool = any(ismember(k, match_ind)) ; 
          
          if catalogFlag
             analysisLogStruct(cc).driver = driver ;  
             analysisLogStruct(cc).effector = effector ;  
             analysisLogStruct(cc).mn_name = mn_name ;  
          end
          cc = cc + 1 ; 
       end
       
    end
    
    % --------------------------------------------------
    %% save all results from one experiment to a sheet
    analysisLogTable = struct2table(analysisLogStruct) ; % convert to table for easy saving
    sheetName = ['Expr' num2str(ExprNum)] ; 
    writetable(analysisLogTable, saveFilename, 'Sheet',sheetName)
    
    fprintf('Completed log for Expr %d \n', ExprNum)
end

% delete extra sheets and turn warning back on
RemoveSheet123(saveFilename)
warning( 'on', 'MATLAB:xlswrite:AddSheet' ) ; % turn off warning