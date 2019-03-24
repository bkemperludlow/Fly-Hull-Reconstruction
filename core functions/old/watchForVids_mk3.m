%Script that watches a directory, listens for creation of .cin files, and
%runs analysis on them as they appear, keeping changelogs of the additions
%mk2: Added queue to replace flags indicating analysis should be done
%mk3: Making sure only the event count stuff is modified by the listener
%events. all other things will take place externally

function watchForVids_mk3()
%pathToWatch = 'D:\Raymond Analysis Test\test automation\' ;
pathToWatch = 'D:\Sam Analysis Test\01_27092016\' ;
ExprNum = 1 ; 
movieCounter = 12 ; %keeps track of movie number in case cameras mess up (number of cine about to save )
pathStruct = generatePathStruct(pathToWatch) ;
cd(pathToWatch)
%3 watchers that look for xy, yz, and xz videos, respectively
fileObjXY = System.IO.FileSystemWatcher(pathToWatch) ;
fileObjXY.EnableRaisingEvents = true ;
fileObjXY.Filter = 'xy*.cin' ;
addlistener(fileObjXY, 'Created', @onChange) ;

fileObjYZ = System.IO.FileSystemWatcher(pathToWatch) ;
fileObjYZ.EnableRaisingEvents = true ;
fileObjYZ.Filter = 'yz*.cin' ;
addlistener(fileObjYZ, 'Created', @onChange) ;

fileObjXZ = System.IO.FileSystemWatcher(pathToWatch) ;
fileObjXZ.EnableRaisingEvents = true ;
fileObjXZ.Filter = 'xz*.cin' ;
addlistener(fileObjXZ, 'Created', @onChange) ;

currEventCount = 0 ; %when all 3 listeners activate analysis will run
eventCount = 0 ; 
currEventCams = cell(1,1) ; %keeps track of which camera's movies are noticed
eventCams = cell(1,1) ; 
currEventNums = [] ; %list of numbers of last 3 videos captured; makes sure numbers agree
eventNums = [] ; 
queue = [] ; %queue for video numbers that need analyzing

while true
    pause(1) 
    if currEventCount >= 3 && mod(currEventCount,3) == 0 %adds to queue once 3 videos (3 cameras) are captured
        add2queue() ; 
    end
    if numel(queue) ~= 0 %do nothing if queue is empty
        %pause(250) %wait for video to be completely saved. No longer
        %needed due to cinFileSize
        runManyAnalysesAuto_mk2(queue(1), ExprNum, pathStruct) 
    end
    queue = queue(2:numel(queue)) ; %remove first element
end
  
    function onChange(~, evt)
        %writes name of file and time registered into the appropriate
        %changelog
        eventName = char(evt.Name) ;
        change = strcat(eventName, ': ', datestr(datetime('now'))) ;
        fileID = fopen(strcat(pathToWatch,'\',eventName(1:2),'log.txt'),'a+') ;
        fprintf(fileID, '%s\r\n', change) ;
        fclose(fileID) ;
        if currEventCount < 3
            currEventCount = currEventCount + 1 ;
            currEventNums = [currEventNums str2num(eventName(4:6))] ;
            currEventCams{currEventCount} = eventName(1:2) ;       
        else
            eventCount = eventCount + 1 ; 
            eventNums = [eventNums str2num(eventName(4:6))] ;
            eventCams{eventCount} = eventName(1:2) ;   
        end
        %disp(eventCount)
        
    end
    
    function add2queue()
        
        %combine xy, xz, and yz logs into one changelog and delete them
        %mergelogs
        %xydir = dir('xy*.cin') ; xzdir = dir('xz*.cin') ; yzdir = dir('yz*.cin') ;
        movieCounter = movieCounter + 1 ;
        camCheck = cellfun(@strcmp, sort(currEventCams(1:3)), {'xy', 'xz', 'yz'}) ;
        
        %makes sure cin files are completely saved before adding to
        %queue
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
        
        %checks that numbers and cameras all match up with expected
        if sum(range(currEventNums)) ~= 0 || ~isequal(camCheck, [1 1 1]) || currEventNums(1) ~= movieCounter
            %if there is a discrepency, change the latest videos and xml files to
            %match movieCounter
            disp('Error: Check video numbers and/or cameras')
            disp(['Using movie number: ' num2str(movieCounter)])
            xydir2 = dir('xy*.xml') ; xzdir2 = dir('xz*.xml') ; yzdir2 = dir('yz*.xml') ;
            movNumStr = fileNumAsStr(movieCounter) ;
            dos(['ren ' xydir(length(xydir)).name ' xy_' movNumStr '.cin']) ;
            dos(['ren ' xzdir(length(xzdir)).name ' xz_' movNumStr '.cin']) ;
            dos(['ren ' yzdir(length(yzdir)).name ' yz_' movNumStr '.cin']) ;
            dos(['ren ' xydir2(length(xydir2)).name ' xy_' movNumStr '.xml']) ;
            dos(['ren ' xzdir2(length(xzdir2)).name ' xz_' movNumStr '.xml']) ;
            dos(['ren ' yzdir2(length(yzdir2)).name ' yz_' movNumStr '.xml']) ;
        end
        %disp('registered 3 events')
        %last_num = eventNums(1) ;
        
        %add video number to queue and update eventCount
        currEventCount = min(eventCount, 3) ; 
        try 
            currEventNums = eventNums(1:min(eventCount, 3)) ;
            currEventCams = eventCams(1:min(eventCount, 3)) ;
        catch
            currEventNums = [] ; 
            eventCams = cell(1,1) ;
        end
        
        eventNums = eventNums(min(eventCount, 3)+1:end) ; 
        eventCams = eventCams(min(eventCount, 3)+1:end) ;
        eventCount = eventCount - min(eventCount, 3) ; 
        
        queue = [queue movieCounter] ;
        disp('Queue:')
        disp(queue)
        %eventCount = eventCount - 3 ;
        %eventCams = eventCams(4:end) ;
        
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