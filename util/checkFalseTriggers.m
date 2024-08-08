% -------------------------------------------------------------------------
% quick function for checking whether or not movies are false triggers --
% useful in cases where we suspect many false triggers are occuring, since
% the normal analysis script takes a while to complete this process
%
% *** update -- for the moment, i can't get cine files to close properly, 
%           so using a hacky workaround
%
% EXAMPLE USAGE:
%{
rootPath = 'D:\Box Sync Old\VNC MN Chrimson\84_24102021\' ;
MovNums = 174:438 ; 

checkFalseTriggers(rootPath, MovNums) ;
%}
% -------------------------------------------------------------------------
function falseTrigStruct = checkFalseTriggers(rootPath, MovNums)
%% --------------------
% inputs + params
if ~exist('MovNums','var') || isempty(MovNums)
    cinDir = dir(fullfile(rootPath, '*.cin*')) ;
    movNumExpression = '(?<camName>[xyz]{2})_(?<movieNum>\d+)' ;
    movNumStrs = arrayfun(@(x) regexp(x.name, movNumExpression,'names'), ...
        cinDir, 'UniformOutput',0) ;
    movieNumbers = cellfun(@(y) str2double(y.movieNum), movNumStrs) ;
    MovNums = unique(movieNumbers) ;
end

% generate path structure for current movie
pathStruct = generatePathStruct(rootPath) ;

% window of time to look into
tWindow = -25:25 ;

% initialize output *** for hacky fix
falseTrigStruct = struct ; 
cc = 1 ; 

%% -------------------------------
% loop over movie numbers
for m = 1:length(MovNums)
    % current movie number
    movNum = MovNums(m) ;
    
    % get directory info for current cine files
    cinDir =  dir(fullfile(rootPath, sprintf('*_%03d.cin*',movNum))) ;
    
    % make sure we have a movie triplet
    if isempty(cinDir)
        fprintf('No files for movie %03d -- skipping ... \n', movNum)
        continue
    elseif length(cinDir) ~= 3
        fprintf('Missing cine file(s) for movie %03d -- skipping ... \n', ...
            MovNum)
        moveToFalseTrig(movNum, pathStruct)
        continue
    end
    
    %% ----------------------------------------------------------
    % get "fly" pixel count in tWindow for each movie (in parallel?)
    pixCountArray = zeros(length(tWindow), length(cinDir)) ;
    
    for n = 1:length(cinDir)
        %% ----------------------
        % get cine path/data
        cinePath = fullfile(cinDir(n).folder, cinDir(n).name) ;
        cineData = myOpenCinFile(cinePath) ; 
        
        % get pixel count
        pixCountArray(:,n) = getCheckFrames(cinePath, cineData, tWindow) ;
        
        % try to close cine file
        myCloseCinFile(cineData) ; 
    end
    
    %% ------------------------------------------
    % check if all movies have adequate pixels
    medianPixCount = median(pixCountArray,1) ; 
    
    if any(medianPixCount < 1)
        fprintf('Movie %03d is a potential false trigger -- should move \n',...
            movNum)
        moveToFalseTrig(movNum, pathStruct) ;
        
%         % *** hacky fix -- store info about movies to move later (after
%         % closing matlab instance)
%         cinDir =  dir(fullfile(pathStruct.root, ...
%             sprintf('*_%03d.cin*',movNum))) ;
%         xmlDir =  dir(fullfile(pathStruct.root, ...
%             sprintf('*_%03d.xml',movNum))) ;
%         combDir = [cinDir ; xmlDir] ;
%         fnList = arrayfun(@(x) fullfile(x.folder, x.name), combDir,...
%             'UniformOutput',false) ; 
%         
%         falseTrigStruct(cc).MovNum = movNum ; 
%         falseTrigStruct(cc).fnList = fnList ; 
%         
%         % increment counter
%         cc = cc + 1 ; 
    else
        fprintf('Movie %03d seems okay -- keeping \n', movNum)
    end
end

% save hacky bullshit struct
save(fullfile(rootPath, 'falseTrigStruct.mat'), 'falseTrigStruct') 

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTION(S)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% -----------------------------------------------------------------------
% 1) get movie bg and 2) test frames around t=0 to look for a fly
function pixCounts = getCheckFrames(cinePath, cineData, tWindow)
%% ------------------
% params
bgNSample = 100 ;
bgMethod = 'Single movie max better';
bwAreaFiltSize = 50 ;

%% -------------------
% get backgound
bg = CineToSparseFormat.FindCineBG(cinePath,bgNSample,bgMethod);

%% ------------------------------------------------------------------
% initialize output (detected "fly" pixels in each frame in tWindow)
pixCounts = zeros(size(tWindow)) ;

%% ---------------------------------------------------------
% loop over tWindow; use method from compression function to find fly pix
for ind = 1:length(tWindow)
    % read current image
    t = tWindow(ind) ;
    inpIm = myReadCinImage(cineData, t);
    
    % apply mask
    mask=imbinarize(bg-inpIm,0.05);
    
    % filter out small objects
    mask=bwareaopen(mask,bwAreaFiltSize); % made for multiple insects in frame
    
    %% --------------------------------------------------------------------
    % fix zeros inside mask to store them in the sparse matrix
    % (black but nonzero!)
    
    inpIm(mask&(inpIm==0))=1;
    
    %% --------------------------------------------------------------------
    % count only what differs from background and store in array
    
    [~,~,v] = find(cast(mask, class(inpIm)).*inpIm);
    pixCounts(ind) = length(v) ; 
    
end

end