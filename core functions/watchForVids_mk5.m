%Script that watches a directory, listens for creation of .cin files, and
%runs analysis on them as they appear, keeping changelogs of the additions
%mk2: Added queue to replace flags indicating analysis should be done
%mk3: Making sure only the event count stuff is modified by the listener
%   events. all other things will take place externally
%mk4: restructuring to have different counts for each cam. should then
%   avoid problems when only two cams save (which happens occasionally)

function watchForVids_mk5()
%pathToWatch = 'D:\Raymond Analysis Test\test automation\' ;
% pathToWatch = 'D:\Sam Analysis Test\01_23062017\';
pathToWatch = 'D:\Box Sync Old\VNC Motor Lines\50_20122018\';

% try to get experiment number automatically
if ~strcmp(pathToWatch(end),'\')
    pathToWatch = [pathToWatch '\'] ; 
end
pathSplit = strsplit(pathToWatch,'\') ; 
folderSplit = strsplit(pathSplit{end-1},'_') ; 
ExprNum = str2double(folderSplit{1}) ; 

mp4Prefix = ['Expr_' num2str(ExprNum)] ; 

clust = parcluster('myLocalProfile') ; 
clust.JobStorageLocation = pathToWatch ; 
clust.HasSharedFilesystem = true ; 
analysis_job_cell = cell(1) ; 
%mp4_job_cell = cell(1) ; 

movieCounter = initializeMovieCounter() ; %keeps track of movie number in case cameras mess up (number of cine about to save )
pathStruct = generatePathStruct(pathToWatch) ;
cd(pathToWatch)

%3 watchers that look for xy, yz, and xz videos, respectively
fileObjXY = System.IO.FileSystemWatcher(pathToWatch) ;
fileObjXY.EnableRaisingEvents = true ;
fileObjXY.Filter = 'xy*.cin' ;
addlistener(fileObjXY, 'Created', @onChangeXY) ;

fileObjYZ = System.IO.FileSystemWatcher(pathToWatch) ;
fileObjYZ.EnableRaisingEvents = true ;
fileObjYZ.Filter = 'yz*.cin' ;
addlistener(fileObjYZ, 'Created', @onChangeYZ) ;

fileObjXZ = System.IO.FileSystemWatcher(pathToWatch) ;
fileObjXZ.EnableRaisingEvents = true ;
fileObjXZ.Filter = 'xz*.cin' ;
addlistener(fileObjXZ, 'Created', @onChangeXZ) ;


xyEventCount = 0 ; xzEventCount = 0 ; yzEventCount = 0 ; 
xyEventNums = [] ; xzEventNums = [] ; yzEventNums = [] ; 
xyEventTimes = [] ; xzEventTimes = [] ; yzEventTimes = [] ; 
%cc = 1 ; 

%currEventCount = 0 ; %when all 3 listeners activate analysis will run
%eventCount = 0 ; 
%currEventCams = cell(1,1) ; %keeps track of which camera's movies are noticed
%eventCams = cell(1,1) ; 
%currEventNums = [] ; %list of numbers of last 3 videos captured; makes sure numbers agree
%eventNums = [] ; 
queue = [] ; %queue for video numbers that need analyzing

while true
    pause(1) 
    if (xyEventCount > 0) && (xzEventCount > 0) && (yzEventCount > 0) %adds to queue once 3 videos (3 cameras) are captured
        add2queue() ; 
    end
    if numel(queue) ~= 0 %do nothing if queue is empty
        %pause(250) %wait for video to be completely saved. No longer
        %needed due to cinFileSize
        %combineAllMoviesInFolder(pathToWatch, queue(1), mp4Prefix) ; 
        %movefile([pathToWatch '*.mp4'], pathStruct.mp4) ; 
        %mp4_job_cell{queue(1)+1} = batch(clust,'combineAllMoviesInFolder',...
        %    0,{pathToWatch, queue(1), mp4Prefix},'CaptureDiary',false) ;
        analysis_job_cell{queue(1)+1} = batch(clust,'flyAnalysisMain',...
            0,{queue(1), ExprNum, pathStruct, true},'CaptureDiary',false) ; 
        %runManyAnalysesAuto_mk2(queue(1), ExprNum, pathStruct) 
    end
    queue = queue(2:numel(queue)) ; %remove first element
    
    completedJobs = findJob(clust,'State','finished') ; 
    if ~isempty(completedJobs)
        %keyboard ;
        delete(completedJobs) ; 
    end
end
  
%----------------------------------------------------------------------------------------------------------------
    function onChangeXY(~, evt)
        %writes name of file and time registered into the appropriate
        %changelog
        eventName = char(evt.Name) ;
        change = strcat(eventName, ': ', datestr(datetime('now'))) ;
        fileID = fopen(strcat(pathToWatch,'\',eventName(1:2),'log.txt'),'a+') ;
        fprintf(fileID, '%s\r\n', change) ;
        fclose(fileID) ;
        
        xyEventCount = xyEventCount + 1 ;
        xyEventNums = [xyEventNums str2num(eventName(4:6))] ;
        xyEventTimes = [xyEventTimes, datetime] ;
    end
    
    function onChangeXZ(~, evt)
        %writes name of file and time registered into the appropriate
        %changelog
        eventName = char(evt.Name) ;
        change = strcat(eventName, ': ', datestr(datetime('now'))) ;
        fileID = fopen(strcat(pathToWatch,'\',eventName(1:2),'log.txt'),'a+') ;
        fprintf(fileID, '%s\r\n', change) ;
        fclose(fileID) ;
        
        xzEventCount = xzEventCount + 1 ;
        xzEventNums = [xzEventNums str2num(eventName(4:6))] ;
        xzEventTimes = [xzEventTimes, datetime] ;
    end

    function onChangeYZ(~, evt)
        %writes name of file and time registered into the appropriate
        %changelog
        eventName = char(evt.Name) ;
        change = strcat(eventName, ': ', datestr(datetime('now'))) ;
        fileID = fopen(strcat(pathToWatch,'\',eventName(1:2),'log.txt'),'a+') ;
        fprintf(fileID, '%s\r\n', change) ;
        fclose(fileID) ;
        
        yzEventCount = yzEventCount + 1 ;
        yzEventNums = [yzEventNums str2num(eventName(4:6))] ;
        yzEventTimes = [yzEventTimes, datetime] ;
    end
%----------------------------------------------------------------------------------------------------------------

    function add2queue()
        
        %combine xy, xz, and yz logs into one changelog and delete them
        %mergelogs
        
        %camCheck = cellfun(@strcmp, sort(currEventCams(1:3)), {'xy', 'xz', 'yz'}) ;
        
        %makes sure cin files are completely saved before adding to
        %queue
        pause(1) ; 
        cinFileSize = 474537120 ;
        xyBytes = 0; xzBytes = 0; yzBytes = 0 ;
        while (xyBytes < cinFileSize || ...
                xzBytes < cinFileSize || ...
                yzBytes < cinFileSize)
            pause(1)
            xydir = dir([pathToWatch 'xy*.cin']) ;
            xzdir = dir([pathToWatch 'xz*.cin']) ;
            yzdir = dir([pathToWatch 'yz*.cin']) ;
            xyBytes = xydir(length(xydir)).bytes ;
            xzBytes = xzdir(length(xzdir)).bytes ;
            yzBytes = yzdir(length(yzdir)).bytes ;
        end
        pause(5)
        
        %checks that movie save times all match up with expected
        currEventTimes = [datenum(xyEventTimes(1)), datenum(xzEventTimes(1)), datenum(yzEventTimes(1))] ;
        eventTimesDist = max(pdist(currEventTimes')) ;
        tol = 5e-4 ; % a little less than a minute
        eventTimesErrorFlag = eventTimesDist > tol ; 
        
        if eventTimesErrorFlag
            disp('Error: mismatch between save times')
            %tol = 1e-4 ; 
            datenumXY = datenum(xyEventTimes) ; datenumXZ = datenum(xzEventTimes) ; datenumYZ = datenum(yzEventTimes) ;   
            coincidenceTimes = intersectTol(datenumXY,intersectTol(datenumXZ,datenumYZ,tol),tol) ; 
            if isempty(coincidenceTimes)
                movieCounter = max([xyEventNums, xzEventNums, yzEventNums]) + 1 ;
                xyEventTimes = [] ; xzEventTimes = [] ; yzEventTimes = [] ;
                xyEventNums = [] ; xzEventNums = [] ; yzEventNums = [] ;
                xyEventCount = 0 ; xzEventCount = 0 ; yzEventCount = 0 ; 
                return ;
            else
                coincidenceTimes = sort(coincidenceTimes) ; 
                xyTimeInd = find(abs(coincidenceTimes(1) - datenumXY) < 2*tol, 1, 'first') ;
                xzTimeInd = find(abs(coincidenceTimes(1) - datenumXZ) < 2*tol, 1, 'first') ;
                yzTimeInd = find(abs(coincidenceTimes(1) - datenumYZ) < 2*tol, 1, 'first') ;
                
                xyEventTimes = xyEventTimes(xyTimeInd:end) ; 
                xzEventTimes = xzEventTimes(xzTimeInd:end) ; 
                yzEventTimes = yzEventTimes(yzTimeInd:end) ;
                
                xyEventNums = xyEventNums(xyTimeInd:end) ; 
                xzEventNums = xzEventNums(xzTimeInd:end) ; 
                yzEventNums = yzEventNums(yzTimeInd:end) ;
                
                xyEventCount = xyEventCount - xyTimeInd + 1 ;
                xzEventCount = xzEventCount - xzTimeInd + 1 ;
                yzEventCount = yzEventCount - yzTimeInd + 1 ;
                
                movieCounter = movieCounter + max([xyTimeInd, xzTimeInd, yzTimeInd]) - 1 ; 
            end
        end
        
        
        %checks that movie numbers all match up with expected
        
        currEventNums = [xyEventNums(1), xzEventNums(1), yzEventNums(1)] ; 
        eventNumErrorFlag = sum(range(currEventNums)) ~= 0  || currEventNums(1) ~= movieCounter ;
        
        if eventNumErrorFlag     
            %if there is a discrepency, change the latest videos and xml files to
            %match movieCounter
            disp('Error: Check video numbers and/or cameras')
            disp(['Using movie number: ' num2str(movieCounter)])
            %{
            xmlFileSize = 204900 ;
            xyBytes = 0; xzBytes = 0; yzBytes = 0 ;
            while (xyBytes < xmlFileSize || ...
                    xzBytes < xmlFileSize || ...
                    yzBytes < xmlFileSize)
                pause(1)
                xydir2 = dir([pathToWatch 'xy*.xml']) ;
                xzdir2 = dir([pathToWatch 'xz*.xml']) ;
                yzdir2 = dir([pathToWatch 'yz*.xml']) ;
                xyBytes = xydir(length(xydir2)).bytes ;
                xzBytes = xzdir(length(xzdir2)).bytes ;
                yzBytes = yzdir(length(yzdir2)).bytes ;
            end
            pause(5) 
            %}
            pause(5)
            
            xydir2 = dir([pathToWatch 'xy*.xml']) ; 
            xzdir2 = dir([pathToWatch 'xz*.xml']) ; 
            yzdir2 = dir([pathToWatch 'yz*.xml']) ;
            
            movNumStr = fileNumAsStr(movieCounter) ;
            dos(['ren "' pathToWatch xydir(length(xydir)).name '" "xy_' movNumStr '.cin"']) ;
            dos(['ren "' pathToWatch xzdir(length(xzdir)).name '" "xz_' movNumStr '.cin"']) ;
            dos(['ren "' pathToWatch yzdir(length(yzdir)).name '" "yz_' movNumStr '.cin"']) ;
            dos(['ren "' pathToWatch xydir2(length(xydir2)).name '" "xy_' movNumStr '.xml"']) ;
            dos(['ren "' pathToWatch xzdir2(length(xzdir2)).name '" "xz_' movNumStr '.xml"']) ;
            dos(['ren "' pathToWatch yzdir2(length(yzdir2)).name '" "yz_' movNumStr '.xml"']) ;
        end
        %disp('registered 3 events')
        %last_num = eventNums(1) ;
        
        %add video number to queue and update eventCounts, eventNums, and
        %eventTimes
        xyEventTimes = xyEventTimes(2:end) ; xzEventTimes = xzEventTimes(2:end) ; yzEventTimes = yzEventTimes(2:end) ;
        xyEventNums = xyEventNums(2:end) ; xzEventNums = xzEventNums(2:end) ; yzEventNums = yzEventNums(2:end) ;
        xyEventCount = xyEventCount - 1 ; xzEventCount = xzEventCount - 1 ; yzEventCount = yzEventCount - 1 ;
        
        queue = [queue movieCounter] ;
        disp('Queue:')
        disp(queue)
        
        movieCounter = movieCounter + 1 ;
       
    end

    %change video number into appropriate string
    function str = fileNumAsStr(fileNum)
        if fileNum < 10
            str = ['00' num2str(fileNum)] ;
        elseif fileNum < 100
            str = ['0' num2str(fileNum)] ;
        else
            str = num2str(fileNum) ;
        end
    end
    
    %get movie number from movie filename
    function fileNum = fileStrAsNum(fileStr)
       if iscell(fileStr)
           fileStr = fileStr{1} ; 
       end
       fileStr_split = strsplit(fileStr,{'_','.'}) ;
       numInd = arrayfun(@(x) ~isnan(str2double(x)), fileStr_split) ; 
       fileNum = str2double(fileStr_split{numInd}) ; 
    end

    function movCount = initializeMovieCounter()
        datadir = dir([pathToWatch '*.cin']) ;
        movNames = {datadir(:).name} ;
        if isempty(movNames)
            movCount = 0 ;
        else
            movNums = arrayfun(@(x) fileStrAsNum(x),movNames) ;
            movCount = max(movNums) + 1 ; 
        end
    end

    function mergelogs
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
end