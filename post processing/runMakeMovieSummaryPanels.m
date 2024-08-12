% -------------------------------------------------------------------------
% script to run "getMovieSummaryPanel.m" on multiple experiment folders
%
% will make this better later
% -------------------------------------------------------------------------
% path info
rootPath = 'D:\Box Sync Old\Test\' ; %'D:\Fly Data\Test\' ;
exprFolders = {} ; % either give full experiment folder name (XX_DDMMYYYY) or just number
exprNumbers = 29 ; 
searchStr = 'Expr_(?<exprNum>\d+)_mov_(?<movieNum>\d+)' ;

% which analysis directories to look in
dirFieldNames = {'unsorted', 'pitchUp', 'pitchDown', 'rollLeft',...
    'rollRight'} ;

% make new summary panel for a video?
overWriteFlag = false ;

% ---------------------------------------------------------------------
% if just experiment numbers are provided, get experiment folder names
if isempty(exprFolders) && ~isempty(exprNumbers)
   % get directory in rootPath. take only folders and remove '..' and '.'
   exprDir = dir(rootPath) ; 
   exprDir = exprDir([exprDir(:).isdir]) ; 
   exprDir = exprDir(3:end) ; 
   
   % get numbers for each experiment folder in dir
   folderStr = '(?<exprNum>\d{2,3})_(?<dateNum>\d+)' ; 
   outNamesAll = arrayfun(@(x) regexp(x.name, folderStr, 'names'), ...
       exprDir, 'UniformOutput', 0) ;
   exprNumsAll = cellfun(@(y) str2double(y.exprNum), outNamesAll) ; 
   
   % find matches between exprNumbers (input) and all experiment numbers;
   % put these into cell
   exprFolders = cell(length(exprNumbers),1) ; 
   for m = 1:length(exprNumbers)
       match_ind = find(exprNumsAll == exprNumbers(m)) ; 
       if ~isempty(match_ind)
            exprFolders{m} = exprDir(match_ind).name ;
       else
           exprFolders{m} = [] ; 
       end
   end
   
elseif isempty(exprFolders) && isempty(exprNumbers)
    fprintf('Need to supply either experiment folder names or numbers \n')
    keyboard
    
end

% -----------------------------------------
% loop over experiment folders
for k = 1:length(exprFolders)
    % check that current exprFolder is not empty
    if isempty(exprFolders{k})
       continue ;  
    end
    
    % current experiment folder
    dataPath = fullfile(rootPath,  exprFolders{k}) ;
    
    % get path struct
    pathStruct = generatePathStruct(dataPath) ;
    
    % ---------------------------------------------------
    % get info for current experiment
    exprInfoStruct = getExprInfo(pathStruct) ;
    cineRangeFrames = [exprInfoStruct.cineFrameStart, ...
        exprInfoStruct.cineFrameEnd] ;
    
    pulseTimingMag = [exprInfoStruct.magPulseStart, ...
        exprInfoStruct.magPulseEnd] ;
    pulseTimingOpto = [exprInfoStruct.optoPulseStart, ...
        exprInfoStruct.optoPulseEnd] ;
    pulseTiming = [pulseTimingOpto ; pulseTimingMag] ;
    
    % ------------------------------------------------------
    % establish where summary panels will be saved
    savePath = fullfile(dataPath, 'summary_panels') ;
    if ~exist(savePath,'dir')
        mkdir(savePath)
    end
    
    % -------------------------------------------
    % loop over analysis folders to look for movie
    for m = 1:length(dirFieldNames)
        % get current data directory
        dataDir = dir(pathStruct.(dirFieldNames{m})) ;
        dataDir = dataDir([dataDir(:).isdir]) ; % take only dir entries
        dataDir = dataDir(3:end) ; % skip '.' and '..'
        
        % loop over movies in current data directory
        for n = 1:length(dataDir)
            tic
            % current movie folder
            folderCurr = dataDir(n).name ;
            
            % save name for summary panel
            savePathFull = fullfile(savePath, ...
                sprintf('%s_summary.png',folderCurr)) ;
            if exist(savePathFull,'file') && ~overWriteFlag
                fprintf('Already created summary for %s \n', folderCurr)
                continue
            end
            
            % get movie number
            out_names = regexp(folderCurr, searchStr,'names') ;
            try
                MovNum = str2double(out_names.movieNum);
            catch
                fprintf('Error getting movie number for %s \n', folderCurr)
                continue
            end
            
            % make summary panel for current movie
            [h_sum, ~] = makeMovieSummaryPanel(pathStruct, MovNum,...
                cineRangeFrames, pulseTiming) ;
            
            % make sure it worked
            if isempty(h_sum)
                continue 
            end
            
            % save summary panel
            try
                exportgraphics(h_sum, savePathFull, 'Resolution', 300) ;
            catch
                print(h_sum, savePathFull, '-dpng','-r300')
            end
            
            % close figure window
            close(h_sum)
            
            toc
            fprintf('Completed %s (%d/%d) \n', folderCurr, n, ...
                length(dataDir))
        end
    end
    
end