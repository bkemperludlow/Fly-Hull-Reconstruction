%Script that watches a directory, listens for creation of .cin files, and
%runs analysis on them as they appear, keeping changelogs of the additions
%mk2: Added queue to replace flags indicating analysis should be done
%mk3: Making sure only the event count stuff is modified by the listener
%   events. all other things will take place externally
%mk4: restructuring to have different counts for each cam. should then
%   avoid problems when only two cams save (which happens occasionally)

function watchForVids_mk9(pathToWatch, cinFileSize)
% ------------------------------------------------------
% check input(s)
if ~exist('pathToWatch','var') || isempty(pathToWatch)
    pathToWatch = 'D:\Box Sync Old\Opto Mechanical\96_04032023\' ; %pwd ;
end

% search expression for movie file names
fnExp = ['(?<camName>[xyz]{2})_(?<movieNum>\d+)|' ...
    '(?<camName>[xyz]{2})_(?<dayNameStr>\w+) (?<monthStr>\w+) ' ...
    '(?<day>\d+) (?<year>\d+) (?<hour>\d+) (?<minute>\d+) ' ...
    '(?<second>\d+.\d+)|'...
    '(?<camName>[xyz]{2})_Y(?<year>\d{4})(?<month>\d{2})(?<day>\d{2})'...
    'H(?<hour>\d{2})(?<minute>\d{2})(?<second>\d+.\d+)']; 

% -------------------------------------------------
% initialize global var values
xyEventCount = 0 ; xzEventCount = 0 ; yzEventCount = 0 ;
xyEventTimes = [] ; xzEventTimes = [] ; yzEventTimes = [] ;
xyFileNames = {} ;  xzFileNames = {} ;  yzFileNames = {} ;

queue = [] ; %queue for video numbers that need analyzing
movFileExt = [] ;

% -------------------------------------------------------------
% try to determine params for checking cine save status
% (this will determine how long we wait for files to save)

% estimate cine file size (either based on calibration or alternative)
cinFileSizeDefault = 1.0e9 ; % hard-coded estimate
if ~exist('cinFileSize','var') || isempty(cinFileSize)
    if exist(fullfile(pathToWatch,'calibration'),'dir')
        calibDir = dir(fullfile(pathToWatch,'calibration','*.cin*')) ;
        if ~isempty(calibDir)
            cinFileSize = max([calibDir.bytes]) ;
        else
            cinFileSize = cinFileSizeDefault ;
        end
    else
        cinFileSize = cinFileSizeDefault ; % expected file size (should make this an input variable up top or something...)
    end
end

% also make an estimate on save speed -- this was just determined roughly,
% and will change based on eg speed of ethernet connection
saveSpeed = 1.0e8 ; % bytes per 5 seconds (5 seconds is the pause in while loop)
byteCountMax = round(2.5*(cinFileSize/saveSpeed)) ;  % how many loops will we allow?

% -----------------------------------
% define cluster
try
    clust = parcluster('myLocalCluster') ;
catch
    clust = parcluster('local') ;
end
clust.JobStorageLocation = pathToWatch ;
clust.HasSharedFilesystem = true ;
analysis_job_cell = cell(1) ;
%mp4_job_cell = cell(1) ;

% get save path info
pathStruct = generatePathStruct(pathToWatch) ;
ExprNum = pathStruct.ExprNum ;
cd(pathToWatch)

% --------------------------------------------------------------------
% define 3 watchers that look for xy, yz, and xz videos, respectively
fileObjXY = System.IO.FileSystemWatcher(pathToWatch) ;
fileObjXY.EnableRaisingEvents = true ;
fileObjXY.Filter = 'xy*.cin*' ;  % ['xy*' movFileExt] ;
addlistener(fileObjXY, 'Created', @(src, evt) onChange(src, evt)) ;

fileObjXZ = System.IO.FileSystemWatcher(pathToWatch) ;
fileObjXZ.EnableRaisingEvents = true ;
fileObjXZ.Filter = 'xz*.cin*' ; % ['xz*' movFileExt] ;
addlistener(fileObjXZ, 'Created', @(src, evt) onChange(src, evt)) ;

fileObjYZ = System.IO.FileSystemWatcher(pathToWatch) ;
fileObjYZ.EnableRaisingEvents = true ;
fileObjYZ.Filter = 'yz*.cin*' ; % ['yz*' movFileExt];
addlistener(fileObjYZ, 'Created', @(src, evt) onChange(src, evt)) ;

fprintf(['\n ' repmat('%',1,100) '\n' ])
fprintf('\n Running experiment: \n %s \n', pathToWatch)
fprintf(['\n ' repmat('%',1,100) '\n'])
% -------------------------------------------------------------------
% MAIN SEQUENCE
% -------------------------------------------------------------------
while true
    % brief pause
    pause(1)

    % -----------------------------------------------------------------
    % if each file listener has seen a file be created, add this
    % triplet to queue
    if (xyEventCount > 0) && (xzEventCount > 0) && (yzEventCount > 0)
        % give a little time (after files are created) to allow all movies
        % to be saved
        pause(10)
        
        % get filenames for matching triplet of movies (also check that
        % there are three movies with matching times)
        [proceedFlag, xy_fn, xz_fn, yz_fn] = groupMovieTriplet(xyEventTimes,...
            xzEventTimes, yzEventTimes, xyFileNames, xzFileNames, ...
            yzFileNames) ;

        %adds to queue once 3 videos (3 cameras) are captured
        if proceedFlag
            add2queue(xy_fn, xz_fn, yz_fn) ;
        else
            % do nothing
            fprintf('No matching movie triplet to add to queue \n')
        end
    end

    % --------------------------------------------------------------
    % if we have something in the queue, run analysis
    if numel(queue) ~= 0 %do nothing if queue is empty
        % if we have something in queue, run analysis
        analysis_job_cell{queue(1)+1} = batch(clust,'flyAnalysisMain',...
            0,{queue(1), ExprNum, pathStruct, true},'CaptureDiary',false,...
            'Pool', [1,4]) ;
        
        % tag analysis job with movie number
        analysis_job_cell{queue(1)+1}.Tag = num2str(queue(1),'%03d') ; 
        
        % after running analysis on most recent queue item (sent to batch),
        % remove from queue
        queue = queue(2:numel(queue)) ; %remove first element
    end
    
    % ----------------------------------------------------------------
    % delete completed jobs to free up cluster resources
    completedJobs = findJob(clust,'State','finished') ;
    if ~isempty(completedJobs)
        fprintf('[DONE] Completed analysis for movie %s (%s)\n', ...
            completedJobs.Tag, datetime(now,'ConvertFrom','datenum'))
        delete(completedJobs) ;
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CALLBACK FUNCTIONS FOR FILE LISTENERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function onChange(~, evt)
    % -----------------------------------------------------
    % gets input from current file being saved, increments
    % movie-specific counters, and writes to log
    % -----------------------------------------------------
    % search filename for movie info
    eventName = char(evt.Name) ;
    out_names = regexp(eventName, fnExp,'names') ;
    
    % if we haven't already specified cine file extension, get it here
    if ~exist('movFileExt','var') || isempty(movFileExt)
        [~, ~, movFileExt] = fileparts(eventName) ; 
    end
    % read out camera name (should be present regardless of file
    % format)
    camName = out_names.camName ;

    % -------------------------------------------------------
    % get datenum information
    if isempty(out_names.movieNum) && isempty(out_names.monthStr)
        % if file includes trigger time (in PCC v3.6 version) use that
        trigger_datestr = strjoin({out_names.year, out_names.month, ...
            out_names.day, out_names.hour, out_names.minute,...
            out_names.second(1:6)},' ') ;
        trigger_datenum = datenum(trigger_datestr, ...
            'yyyy mm dd HH MM SS.FFF') ;
    elseif isempty(out_names.movieNum) && ~isempty(out_names.monthStr)
        % if file format includes trigger time (in "~5" format), use that
        trigger_datestr = strjoin({out_names.year, out_names.monthStr, ...
            out_names.day, out_names.hour, out_names.minute,...
            out_names.second},' ') ;
        trigger_datenum = datenum(trigger_datestr, ...
            'yyyy mmm dd HH MM SS.FFF') ;
    else
        % otherwise, use current datetime as a proxy for trigger time
        % NB: THIS IS LAZY -- SHOULD GET FROM XML FILE
        trigger_datenum = datenum(datetime('now')) ;
    end

    % ----------------------------------
    % update log file
    change = strcat(camName, ': ', datestr(trigger_datenum)) ;
    fileID = fopen(strcat(pathToWatch,'\',eventName(1:2),'log.txt'),'a+') ;
    fprintf(fileID, '%s\r\n', change) ;
    fclose(fileID) ;

    % ---------------------------------
    % update global counters
    switch camName
        case 'xy'
            xyFileNames{end+1} = eventName ;
            xyEventTimes = [xyEventTimes; trigger_datenum] ;
            xyEventCount = xyEventCount + 1 ;
        case 'xz'
            xzFileNames{end+1} = eventName ;
            xzEventTimes = [xzEventTimes; trigger_datenum] ;
            xzEventCount = xzEventCount + 1 ;
        case 'yz'
            yzFileNames{end+1} = eventName ;
            yzEventTimes = [yzEventTimes; trigger_datenum] ;
            yzEventCount = yzEventCount + 1 ;
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------------------
% take a set of three movies recently saved, make sure their
% timing matches up, and then output filenames for these movies
%
% NB: taking event times and filenames as inputs to reduce risk of
% list being updated during computation
% ---------------------------------------------------------------------
function [proceedFlag, xy_fn, xz_fn, yz_fn] = groupMovieTriplet(xyTimes,...
        xzTimes, yzTimes, xyNames, xzNames, yzNames)
    
    % define datenum tolerance for matching triplet
    tol = 5e-5 ; % *** WAS 5e-4 (a little more than a minute) but now that camera times are synced i think we can go lower?

    % find all combinations of possible datenums in pre-queue
    [ii, jj, kk] = ndgrid(xyTimes, xzTimes, yzTimes) ;
    eventTimeCombs = [ii(:), jj(:), kk(:)] ;

    % get combination with smallest difference in event timing
    timeDist = max(abs([(eventTimeCombs(:,1) - eventTimeCombs(:,2)), ...
        (eventTimeCombs(:,1) - eventTimeCombs(:,3)), ...
        (eventTimeCombs(:,2) - eventTimeCombs(:,3))]),[],2) ;
    [minTimeDist, minInd] = min(timeDist) ;

    % make sure the total distance between time points is within
    % tolerance
    if (minTimeDist < tol)
        % in this case, distance is within tolerance and we should
        % select this combination
        proceedFlag = true ;

        % find the indices in camera event times array corresponding to
        % the datenum combination with minimal time gap
        [~, xy_ind] = min(abs(xyTimes - eventTimeCombs(minInd,1))) ;
        [~, xz_ind] = min(abs(xzTimes - eventTimeCombs(minInd,2))) ;
        [~, yz_ind] = min(abs(yzTimes - eventTimeCombs(minInd,3))) ;

        % use these indices to get the matching filenames
        xy_fn = xyNames{xy_ind} ;
        xz_fn = xzNames{xz_ind} ;
        yz_fn = yzNames{yz_ind} ;
    else
        proceedFlag = false ;
        xy_fn = [] ; xz_fn = [] ; yz_fn = [] ;
        decrementCounters()
    end

end

% ------------------------------------------------------------
% take a triplet of movie files and add them to the queue for
% analysis
%
% To add: merge logs
% -------------------------------------------------------------
function add2queue(xy_fn, xz_fn, yz_fn)
    % -------------------------------------------------------------
    % make sure cin files are completely saved before adding to
    % queue (i.e. make sure file size is sufficiently large)
    pause(1) ;
    xyBytes = 0; xzBytes = 0; yzBytes = 0 ;
    byteCounter = 0 ; 
    while (xyBytes < cinFileSize || ...
            xzBytes < cinFileSize || ...
            yzBytes < cinFileSize) && (byteCounter < byteCountMax)
        % give a little time to allow save to progress
        pause(5)

        % get current directory for camera movie files
        xydir = dir(fullfile(pathToWatch, ['xy*' movFileExt])) ;
        xzdir = dir(fullfile(pathToWatch, ['xz*' movFileExt])) ;
        yzdir = dir(fullfile(pathToWatch, ['yz*' movFileExt])) ;

        % get index in directory for current movie triplet filenames
        xy_ind = arrayfun(@(x) strcmp(x.name, xy_fn), xydir) ;
        xz_ind = arrayfun(@(x) strcmp(x.name, xz_fn), xzdir) ;
        yz_ind = arrayfun(@(x) strcmp(x.name, yz_fn), yzdir) ;
                
        % check current file size
        xyBytes = xydir(xy_ind).bytes ;
        xzBytes = xzdir(xz_ind).bytes ;
        yzBytes = yzdir(yz_ind).bytes ;
        
        % increment counter to determine how many loops we've gone through
        byteCounter = byteCounter + 1 ; 
    end
    
    % if we broke while loop because of too many loops, it might mean file
    % size doesn't match -- just skip over these for now and deal with them
    % later
    if byteCounter == byteCountMax
        divertData({xy_fn, xz_fn, yz_fn}, {movFileExt, '.xml'})
        return
    end
    
    % ... otherwise pause to allow any final saving processes to take place
    pause(5)

    % -----------------------------------------------------------------
    % if applicable, rename movies from datetime filename format to 3
    % digit integer format
    regexp_check = regexp({xy_fn, xz_fn, yz_fn}, fnExp, 'names') ;

    % first check if (unaltered) files have associated movie numbers
    % already. NB: numberMismatchFlag should correspond to cases wherein 
    % the movie numbers don't match up; assignNumberFlag should correspond
    % to cases wherein the files are being saved with the trigger time in 
    % the filename
    movNums = cellfun(@(y) str2double(y.movieNum), regexp_check) ;
    numberMismatchFlag = any(isnan(movNums)) | ...
        ~all(movNums(:) == movNums(1)) ;
    assignNumberFlag = all(isnan(movNums)) ;
    
    % regardless of reason why, define flag for any renaming that needs to
    % occur
    renameFlag = assignNumberFlag | numberMismatchFlag ;
    
    % generate lists that we can loop over to make process of
    % renaming easier (if we're renaming)
    if renameFlag
        old_fn_list = {xy_fn, xz_fn, yz_fn} ;
        cam_names = {'xy', 'xz', 'yz'} ;
        file_type_list = {movFileExt, '.xml'} ;
    end
    
    % if either we have non-matching movie numbers OR datetime format
    % movie filenames, rename the cine and xml files
    if assignNumberFlag
        % if we need to rename movies, want to get current maximum
        % movie number (to avoid overwriting)
        movieCountCurr = getMovieCounter() ; % current max movie number
        
        % TESTING -- with datetime file format, should have less risk of
        % overwriting, so just use (movieCountCurr + 1) as new movie number
        movieCounter = movieCountCurr + 1 ;
        
    elseif numberMismatchFlag
        % if we need to rename movies, want to get current maximum
        % movie number (to avoid overwriting)
        movieCountCurr = getMovieCounter() ; % current max movie number
        
        % define a range of possible movie numbers we could use, in
        % ascending order
        minMovNum = min(movNums,[],'omitnan') ;
        maxMovNum = movieCountCurr ;
        movNumRange = (minMovNum:maxMovNum) ;
        
        % loop over possible movie numbers and see if we can use any
        % without overwriting data
        movNumCheck = false(size(movNumRange)) ;
        for q = 1:length(movNumRange)
            % for a new movie number to be a good one, we should be
            % able to save all movies to this number without
            % overwriting. so we check that
            nonMatchingCams = cam_names(movNums ~= movNumRange(q)) ;
            movNumStr = num2str(movNumRange(q), '%03d') ;
            movNumCheck(q) = any(cellfun(@(y) ~exist(fullfile(pathToWatch, ...
                [y '_' movNumStr movFileExt]),'file'), nonMatchingCams)) ;
        end
        
        % check if we found a good movie number for all cameras. if so,
        % use that. otherwise move data elsewhere to be analyzed later
        if sum(movNumCheck) > 0
            movieCounter = movNumRange(find(movNumCheck,1,'first')) ;
        else
            % here we need a better solution, but for now i'm just
            % going to shuttle these movies to a different folder to
            % prevent data loss
            %{
            t_now = now ;
            tempSaveDir = fullfile(pathToWatch, num2str(t_now)) ;
            mkdir(tempSaveDir)
            
            % loop over cameras and file types to move data
            for mm = 1:length(cam_names)
                % get basename for old filename (for this camera)
                [~, base_name, ~] = fileparts(old_fn_list{mm}) ;
                for nn = 1:length(file_type_list)
                    % current file extension
                    ext_curr = file_type_list{nn} ;
                    fn_wExt = [base_name ext_curr] ;
                    
                    % move file
                    old_fn = fullfile(pathToWatch, fn_wExt) ;
                    new_fn = fullfile(tempSaveDir, fn_wExt) ;
                    movefile(old_fn, new_fn) ;
                end
            end
            
            % after moving files, exit without adding to queue
            decrementCounters(xy_fn, xz_fn, yz_fn)
            %}
            divertData(old_fn_list, file_type_list)
            return
        end
    else
        % in this case, all movie numbers match, so we're good
        movieCounter = movNums(1) ;
    end
    
    
    % -------------------------------------------------------------
    % now that we have movie number sorted out (hopefully), perform
    % renaming if necessary
    if renameFlag
        for mm = 1:length(cam_names)
            % get current camera name
            cam = cam_names{mm} ;
            % get basename for old filename (for this camera)
            [~, base_name, ~] = fileparts(old_fn_list{mm}) ;
            
            % loop over file types
            for nn = 1:length(file_type_list)
                % current file extension
                ext_curr = file_type_list{nn} ;
                
                % rename file
                old_fn = fullfile(pathToWatch, [base_name ext_curr]) ;
                new_fn = fullfile(pathToWatch, ...
                    [cam '_' num2str(movieCounter,'%03d') ext_curr]) ;
                
                if ~strcmp(old_fn, new_fn)
                    status = movefile(old_fn, new_fn) ;
                    if ~status
                        fprintf('Error changing filename from %s to %s \n',...
                            old_fn, new_fn)
                        decrementCounters(xy_fn, xz_fn, yz_fn)
                        return
                    end
                end
            end
        end
    end

    % ----------------------------------------------------------------
    % add current movie to the queue and decrease event counters/remove
    % entries from event time lists
    queue(end+1) = movieCounter ;
    decrementCounters(xy_fn, xz_fn, yz_fn)
    fprintf('[ADD] Movie %03d added to analysis queue (%s) \n', ...
        movieCounter, datetime(now,'ConvertFrom','datenum'))

end

% --------------------------------------------------------------------
% get maximum number of current movie in directory and return that
% number (or zero, if no movies)
function movCount = getMovieCounter()
    datadir = dir(fullfile(pathToWatch, ['*' movFileExt])) ;
    movNames = {datadir(:).name} ;
    if isempty(movNames)
        movCount = 0 ;
    else
%         searchExp = '(?<camName>[xyz]+)_(?<movieNum>\d+)' ;
        out_names = regexp(movNames, fnExp, 'names') ;
        movNums = cellfun(@(x) str2double(x.movieNum), out_names) ;
        movCount = max([movNums, 0],[],'omitnan') ;
        
    end
end

% ---------------------------------------------------------------------
%% decrease event counters for the three cameras
function decrementCounters(xy_fn, xz_fn, yz_fn)
    % decrease event counter by 1 for each movie
    xyEventCount = xyEventCount - 1 ;
    xzEventCount = xzEventCount - 1 ;
    yzEventCount = yzEventCount - 1 ;
    
    % if we have specific files as inputs, remove the event times and file
    % names
    if (nargin > 0)
        % get indices for filenames
        xy_ind = cellfun(@(x) strcmp(x, xy_fn), xyFileNames) ;
        xz_ind = cellfun(@(x) strcmp(x, xz_fn), xzFileNames) ;
        yz_ind = cellfun(@(x) strcmp(x, yz_fn), yzFileNames) ;
        
        % remove filenames and event numbers, decrement event counter
        xyFileNames = xyFileNames(~xy_ind) ;  % xy
        xyEventTimes = xyEventTimes(~xy_ind) ;
        
        xzFileNames = xzFileNames(~xz_ind) ;  % xz
        xzEventTimes = xzEventTimes(~xz_ind) ;
        
        yzFileNames = yzFileNames(~yz_ind) ;  % yz
        yzEventTimes = yzEventTimes(~yz_ind) ;  
    end
end

% ---------------------------------------------------------------------
%% merge txt log files -- currently unused
function mergelogs()
    xyID = fopen([pathToWatch 'xylog.txt']) ;
    xzID = fopen([pathToWatch 'xzlog.txt']) ;
    yzID = fopen([pathToWatch 'yzlog.txt']) ;
    %changeID = fopen([pathToWatch 'changelog.txt'], 'w+') ;
    changeID = fopen([pathToWatch 'changelog.txt'], 'a+') ;
    xyend = '' ; xzend = '' ; yzend = '' ;
    while ~isnumeric(xyend) && ~isnumeric(xzend) && ~isnumeric(yzend)
        xyend = fgets(xyID) ; xzend = fgets(xzID) ; yzend = fgets(yzID) ;
        fprintf(changeID, '%s\r', xyend) ;
        fprintf(changeID, '%s\r', xzend) ;
        fprintf(changeID, '%s\r', yzend) ;
    end

    fclose(xyID) ; fclose(xzID) ; fclose(yzID) ; fclose(changeID) ;
    delete([pathToWatch 'xylog.txt'], [pathToWatch 'xzlog.txt'], [pathToWatch 'yzlog.txt']) ;
end

% ---------------------------------------------------------------------
%% remove data from main directory due to potential save error
function divertData(fn_list, file_type_list)
    t_now = now ;
    tempSaveDir = fullfile(pathToWatch, num2str(t_now)) ;
    mkdir(tempSaveDir)
    
    % loop over cameras and file types to move data
    for mm = 1:length(fn_list)
        % get basename for old filename (for this camera)
        [~, base_name, ~] = fileparts(fn_list{mm}) ;
        for nn = 1:length(file_type_list)
            % current file extension
            ext_curr = file_type_list{nn} ;
            fn_wExt = [base_name ext_curr] ;
            
            % move file
            old_fn = fullfile(pathToWatch, fn_wExt) ;
            new_fn = fullfile(tempSaveDir, fn_wExt) ;
            movefile(old_fn, new_fn) ;
        end
    end
    
    % after moving files, exit without adding to queue
    decrementCounters(fn_list{1}, fn_list{2}, fn_list{3})
end
end