% -------------------------------------------------------------------------
%% Cine to Sparse converter
% Converts cine files to mat files with frames containing moving objects in 
% a sparse format, and a metadata struct
% 
% From Tsevi Beatus Lab, August 2021  
% Adapted by Sam Whitehead (Aug 2021) to:
%   1. Adapt for older file formats
%   2. Remove mpg creation, which can be added back in, but in original
%   version uses MATLAB classes we don't have 
%
% INPUTS:
%   - filenames: cell array of full filenames for movies to compress
%   - bgMethod: method used to extract background from movie
%   - checkBGFlag: boolean variable which, if true, prompts the user to
%       approve of the calculated background
% 
% OUTPUTS:
%   - N/A
%
% EXAMPLE USAGE:
%{
path = 'D:\Box Sync Old\VNC MN Chrimson\84_24102021\'
fn = 'xy_166.cine' ; 
fn_full = fullfile(path, fn) ; 

tic
CineToSparseFormat.Cine2SparseSam(fn_full) ;
toc
%}
% -------------------------------------------------------------------------
function [] = Cine2SparseSam(filenames, bgMethod, checkBGFlag, ...
    defaultFrameRate, overWriteFlag)
% -------------------
%% inputs and params
if ~exist('filenames','var') || isempty(filenames)
    filterSpec = pwd ; % CHANGE FOR CONVENIENCE
    filenames = uipickfiles('FilterSpec',filterSpec,...
        'Type',{'*.ci*', 'All files' },'Prompt','Select Cine\s to convert') ;
    
    if isempty(filenames)
        disp('User selected Cancel; end of session')
        return
    end
end
if ~exist('bgMethod','var') || isempty(bgMethod)
    bgMethod='Single movie max better';
     %{
    'Single movie max better'
    'Single movie max'
    'Single movie mean'
    'Single movie median'
    'Manual input'
    %}
end
if ~exist('checkBGFlag','var') || isempty(checkBGFlag)
    checkBGFlag = false ; 
end
if ~exist('defaultFrameRate','var') || isempty(defaultFrameRate)
    % default frame rate for cine (in case we can't get it from metadata)
    defaultFrameRate = 8000 ; % frames per second
end
if ~exist('overWriteFlag','var') || isempty(overWriteFlag)
    % create a new sparse struct for movie if one already exists?
    overWriteFlag = false ; 
end
% ----------------------------------
% params
% ----------------------------------
cornell = 1 ; % not sure what the purpose of this is yet -- something about 8 vs 16 bit video

% background calculation params
bgNSample=100;
filterMethod='Area filter';
%{
'Largest blob'
'Area filter'
%}

% valid movie file extentions
validExts = {'.cin', '.cine'} ; 

% movie search str 
movSearchStr = '(?<camName>[xyz]{2})_(?<MovNum>\d{3})' ;

% sparse data save string
saveFnStr = '%s_sparse.mat' ; 

% make sure filenames is cell type
if ~iscell(filenames)
    filenames={filenames};
end
% ------------------------------------------------------------------------
% % Load PhantomSDK if needed: (!!take care: PHANTOM SDK LOADING NEEDS TO 
% % BE REVISED GLOBALLY!!)
% % NB: should be loaded by startup file

% if ~exist('getCinMetaData','file')
%     run('PhantomSDK\runMeFirst.m')
% end

% ---------------------------------------------------------------
%% Loop on cine files in the path
fileInd=0;

while fileInd<length(filenames)
%%
% increment file counter
    fileInd=fileInd+1;
    
%% 
% check file format
    cinePath = filenames{fileInd} ; 
    [path, fn, ext] = fileparts(cinePath) ; 
    if ~ismember(ext, validExts)
        fprintf('Removing %s -- not appropriate file type\n',...
            cinePath)
        continue
    end
    
%%
% check if cine has been compressed already
savePathTest = fullfile(path, sprintf(saveFnStr,fn)) ; 
if exist(savePathTest,'file') && ~overWriteFlag
   fprintf('Already compressed %s -- skipping \n', savePathTest)
   continue
end

%% 
% Load the cine file:
    cineMetaData = getCinMetaData(cinePath) ;
    cineData  = myOpenCinFile(cinePath);
    fprintf('Compressing %s ...\n', cinePath);
    
%% 
% Choose the background generation algorithm:
% SCW: adding option to load background from previously processed data
% struct. assumes movie is in format (cam name)_(3 digit movie num). Try
% this first

    [bg, errorFlag] = loadPrevBG(path, fn, movSearchStr) ; 
    
    if ~errorFlag
        metaData.bg = bg ; 
    else
        switch strtok(bgMethod)
            case 'Manual'
                [file,path] = uigetfile('*.*',...
                    ['Choose background image for ',filenames{fileInd}]);
                if isequal(file,0)
                   disp('User selected Cancel; skipping cine');
                   continue
                else
                   metaData.bg = imread(fullfile(path,file));
                end
            case 'Single'
                metaData.bg = CineToSparseFormat.FindCineBG(cinePath,...
                    bgNSample,bgMethod);
        end
    end
    
%% 
% Generate the background image and check with user if its good:
    if checkBGFlag
        checkFig=figure;
        set(gcf,'Visible','on')
        imshow(metaData.bg,[])

    %     answer='Yes';
        answer = questdlg('Is background good?', ...
            'Background Check', ...
            'Yes','No','Yes');
        % answer = input('Is the background okay? (Yes/No)') ;
        delete(checkFig)

        switch answer
            case 'Yes'
            otherwise
                disp('Background is bad; skipping cine');
                continue
        end
    end
%% 
% Assign and initialize variables

    numOfImages = cineMetaData.lastImage-cineMetaData.firstImage+1 ;
    clear frames
    frames(numOfImages,1)=struct; %#ok<SAGROW>
    metaData.startFrame=cineMetaData.firstImage;
    if isfield(cineMetaData,'frameRate')
        % to get cine frame rate, first try to read from metadata (this
        % should work for videos saved using newer software versions)
        metaData.frameRate=cineMetaData.frameRate;
    else
        % otherwise, try to get frame rate from xml file 
        xml_filename = fullfile(path, [fn, '.xml']) ; 
        if exist(xml_filename, 'file')
           % if xml file exists, load that and get frame rate
           xml_info = getInfoXML(xml_filename) ; 
           metaData.frameRate=xml_info.frameRate;
        else
            % otherwise use default (hard-coded) frame rate
            metaData.frameRate = defaultFrameRate ;
        end
    end
    metaData.frameSize=size(metaData.bg);
   
%% Loop on frames in the cine

    saveInd=0;
    for imageInd = 1:numOfImages
        saveInd=saveInd+1;
        if ~mod(imageInd,50)
            disp([num2str(imageInd),'/',num2str(numOfImages)]);
        end
        inpIm = myReadCinImage(cineData, cineMetaData.firstImage+imageInd-1);
        mask=imbinarize(metaData.bg-inpIm,0.05);
        switch filterMethod
            case 'Largest blob'
                mask=bwareafilt(mask,1); % choose only the largest blob
            case 'Area filter'
                mask=bwareaopen(mask,50); % made for multiple insects in frame
        end
        
%% 
% fix zeros inside mask to store them in the sparse matrix (black but nonzero!)

        inpIm(mask&(inpIm==0))=1;
%% 
% copy only what differs from background and save as sparse matrix

        [row,col,v] = find(cast(mask, class(inpIm)).*inpIm);
        if cornell==1
            v = uint16(v);
            frames(saveInd).indIm=uint16([row,col,uint16(v*255)]);
        else
            frames(saveInd).indIm=uint16([row,col,uint16(v)]);
        end

    end
%% 
% Save xml file and add to metaData

     xml_filename = fullfile(path, [fn '.xml']) ; 
     if exist(xml_filename, 'file')
         % if xml file exists, load that and get struct
         xmlStruct = parseXML(xml_filename) ; 
         
         % remove timeblock/exposureblock info that just eats up space
         xmlStruct.chd.TIMEBLOCK.Date = xmlStruct.chd.TIMEBLOCK.Date{1} ;
         xmlStruct.chd.TIMEBLOCK.Time = xmlStruct.chd.TIMEBLOCK.Time{1} ;
         xmlStruct.chd.EXPOSUREBLOCK.Exp = xmlStruct.chd.EXPOSUREBLOCK.Exp{1} ; 
         
         % add to metaData
         metaData.xmlStruct = xmlStruct.chd ; 
     else
         % otherwise ???
         metaData.xmlStruct = [] ; 
     end  
%% 
% Check that file isn't too large (problems with object segmentation)

    framesWhos=whos('frames');
    if framesWhos.bytes>10^9
        warndlg('File is Extremely large, something is probably wrong!')
    end
%% 
% save the frame and metadata variables:
% 
% (v7.3 is important for partial loading of images)
    save(fullfile(path, sprintf(saveFnStr,fn)), 'frames','metaData','-v7.3')
    myCloseCinFile(cineData);
    fprintf('Completed compression of %s \n', cinePath);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------------------------------------------------------------------
%% function to try to get movie background from previously analyzed results
function [bg, errorFlag] = loadPrevBG(path, fn, movSearchStr)
% -----------------------------------
% initialize output
bg = [] ; 
errorFlag = true ; 

% --------------------------------------------
% get movie number from file name
tokens = regexp(fn, movSearchStr,'names') ; 
if isempty(tokens)
    return
end
MovNum = str2double(tokens.MovNum) ; 
camName = tokens.camName ; 

% -----------------------------------------------------------
% try to find path for movie analysis output
pathStruct = generatePathStruct(path) ; 
movAnalysisPath = findMovAnalysisPath(pathStruct, MovNum) ; 
if isempty(movAnalysisPath)
    return
end

% -----------------------------------------------
% load movie analysis output
[~, analysisOutput, ~, ~, ~] = hierarchicalLoadData(pathStruct, MovNum,...
    [],[],true) ;

% get all bg and params
allBG = analysisOutput.allBG ;
params = analysisOutput.params ; 

% get bg for current camera view
bg = squeeze(allBG(params.(upper(camName)),:,:)) ; 
errorFlag = false ; 
end