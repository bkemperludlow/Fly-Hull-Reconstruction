% -------------------------------------------------------------------------
% quick script to get xml data for movies
% -------------------------------------------------------------------------
rootPath = 'D:\Box Sync Old\Opto Silencing\47_24102020\' ;
saveFlag = true ;

% get directory of all xml files in current directory
xmlDir = dir(fullfile(rootPath, '*.xml')) ;

% initialize output struct
xmlStruct = struct() ;

% loop over files and read out xml info
tic
parfor ind = 1:length(xmlDir)
    % current xml filename
    filename = fullfile(xmlDir(ind).folder, xmlDir(ind).name) ;
    
    % read out xml info in struct form
    xmlInfo = parseXML(filename) ;
    
    % add xml info to main struct
    [xmlStruct(ind).path, xmlStruct(ind).fn, ~] = fileparts(filename) ;
    xmlStruct(ind).info = xmlInfo ;
    
    
    % fprintf('Completed %d/%d xml files \n', ind, length(xmlDir)) 
end
toc

% shut down parallel pool
delete(gcp) ; 

% save results?
if saveFlag
    save(fullfile(rootPath, 'xmlStruct.mat'), 'xmlStruct')
end

% -------------------------------------------------------------
% also grab trigger times/ fps
triggerTimes = cell(3,1) ;
fpsList = cell(3,1) ; 
firstImList = cell(3,1) ; 

prefix_list = {'xy', 'xz', 'yz'} ; 
for j = 1:length(prefix_list)
    % get index for current camera
    idx = arrayfun(@(x) contains(x.fn, prefix_list{j}), xmlStruct) ; 
    
    % get trigger times
    trigCurr = ...
        arrayfun(@(x) datenum(x.info(2).Children(1).Children(8).Children(2).Children.Data(1:12),'HH:MM:SS.FFF'),...
        xmlStruct(idx)) ;
    
    triggerTimes{j} = trigCurr ; 
    
    % get frame rate
    frameRate = ...
        arrayfun(@(x) str2double(x.info(2).Children(3).Children(28).Children.Data),...
        xmlStruct(idx)) ;
    
    fpsList{j} = frameRate ;
    
    % get first frame number
     firstIm = ...
        arrayfun(@(x) str2double(x.info(2).Children(1).Children(6).Children.Data),...
        xmlStruct(idx)) ;
    
    firstImList{j} = firstIm ;
end

% xz_idx = arrayfun(@(x) contains(x.fn, 'xz'), xmlStruct) ; 
% yz_idx = arrayfun(@(x) contains(x.fn, 'yz'), xmlStruct) ; 
% movNums = unique(arrayfun(@(x) str2double(x.fn(4:end)), xmlStruct)) ; 
% 
% triggerTimes = 