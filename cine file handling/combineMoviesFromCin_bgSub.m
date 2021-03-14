% -------------------------------------------------------------------------
% function to generate 3-panel mp4 file from cine files with background
% subtracted and a patch to indicate stimulus period
% -------------------------------------------------------------------------
function stat = combineMoviesFromCin_bgSub(metaData, outputFileName, ...
    output_fps, realFrameRate,allBG, pulseDuration, patchColor )
% -----------------------
%% inputs and params
if ~exist('pulseDuration','var') || isempty(pulseDuration)
    pulseDuration = 50 ; % in milliseconds
end
if ~exist('patchColor','var') || isempty(patchColor)
    patchColor = [] ; % color for 'insertShape'. If empty, won't draw patch
    % NB: can be 'blue' | 'green' | 'red' | 'cyan' | 'magenta' | 'black' | 'black' | 'white'
end

% start timer
tic ;

stat = 0 ;

% get camera indices
XZ = 1 ;
XY = 2 ;
YZ = 3 ;

% index of background images for each camera
BG_ind = [2 3 1] ;

% time between frames
dt = 1/output_fps ;

% frame number size for each video chunk
chunksize = 365 ;

% video frame rate
if (~exist('realFrameRate','var'))
    realFrameRate = 8000 ;
end

% filename for video output
mp4Filename = [outputFileName '.mp4'] ;
%conf.videoQuality = 75;

% coordinates for stimulus patch
patchCoords =  [1 1 35 35; 513 1 35 35 ; 1025 1 35 35 ] ;

% ------------------------------------------------------------
%% initialize timer on the bottom left corner
Ny = 15 ;
Nx = 200 ;
htimer = figure ;
set(htimer,'Renderer','zbuffer')
%set(gcf,'Renderer','painters');

I = zeros(Ny,Nx,'uint8') + uint8(255) ;
figure(htimer) ;

imshow(I);
set(gca,'position',[0 0 1 1 ]) ;
set(gca,'units','pixels') ;
set(gca,'position',[0 0 Nx Ny]) ;
set(gcf,'position',[500 700 Nx Ny]) ;
set(gca,'position',[0 0 Nx Ny]) ;
set(htimer,'paperPositionMode','auto') ;
st = 'bla' ;
htext = text( 'units','pixels','position',[1 0],'fontsize',12,'FontWeight','bold',...
    'VerticalAlignment','bottom','string',st) ;

% ---------------------------------------------------------
%% get first image, last image, and total image number
% find the min firstImage and max lastImage of all movies
firsts = [metaData.firstImage] ;
lasts  = [metaData.lastImage] ;
flags  = ([metaData.exists]==1) ;

firsts_ = firsts(flags) ; % only the ones that exist
lasts_  = lasts(flags) ;

% calculate the total numbers of frames to read, to use in the progress bar
Nprogress = sum(lasts_-firsts_) + numel(firsts_) ;

globalFirstImage = min(firsts_) ;
globalLastimage  = max(lasts_) ;
Ntot = globalLastimage - globalFirstImage + 1 ;

clear firsts_ lasts_

%movObj = cell(1,3) ;

% -----------------------------------------------------------
%% check image dimensions across cameras
cams = [XZ XY YZ] ;
availableCams = cams(flags) ;

Hvec = [metaData.height] ;
Wvec = [metaData.width] ;

Hvec = unique(Hvec(flags)) ;
Wvec = unique(Wvec(flags)) ;

if (numel(Hvec)>1)
    disp('Error: Movies have different image heights.') ;
    stat = -1 ;
    return ;
end
if (numel(Wvec)>1)
    disp('Error: Movies have different image widths.') ;
    stat = -1 ;
    return ;
end

H = Hvec ; W = Wvec ; clear Hvec Wvec ;

% ------------------------------------------------------------------
%% initialize waitbar and video writer
Nchunks = ceil(Ntot/chunksize) ;

t = globalFirstImage ; % current "time" counted in frames
currtime = 0 ;

ttl = ['Processing movie ' mp4Filename ] ;
ttl(ttl=='_')=' ' ;

h = waitbar(0,ttl);

waitbarCounter = 0 ;

writerObj=VideoWriter(mp4Filename,'MPEG-4');
open(writerObj);

% ---------------------------------------------------------------
%% loop through frame chunks and write to video
for chunk=1:Nchunks
    clear vid
    
    if (chunk<Nchunks)
        N = chunksize ; % length(xzInfo) ;
    else
        N = Ntot - chunksize*(chunk-1) ;
    end
    
    % define vid - the image stack to be added to the movie
    vid.times = currtime + (0:(N-1))*dt;
    vid.height = H ;
    vid.width  = W * 3;
    fr = zeros(vid.height, vid.width, 3,'uint8') ;
    vid.frames(1:N)  =  struct('cdata',fr);
    
    % -----------------------------------------------------
    %% loop over cameras
    for k=availableCams
        c = 0 ; % counter into vid.frames
        
        % go over images in this cam and chunk
        for it = t:(t+N-1)
            c=c+1 ;
            % if current image is out of range
            if ( it<firsts(k) || it>lasts(k))
                continue ;
            end
            % read image "it" from camera "k"
            % first find the index of the image in the movie
            n = it - firsts(k) + 1 ;
            try
                % try to process image w/ bg subtraction
                im1 = myReadCinImage(metaData(k).cindata, it) ;
                im2 = imsubtract(squeeze(allBG(BG_ind(k),:,:)),im1) ;
                low_in = double(min(im2(:)))/255 ;
                high_in = double(max(im2(:)))/255 ;
                if high_in <= low_in
                    img = imcomplement(im2) ;
                else
                    im3 = imadjust(im2,[low_in,high_in],[0,1]) ;
                    img = imcomplement(im3) ;
                end
                %img = imsubtract(
            catch ERR
                disp(ERR) ;
                keyboard;
            end
            waitbarCounter = waitbarCounter + 1 ;
            % place image inside triple frame
            % when k=1 (xz) use 1:512 -->   1:W
            % when k=2 (xy) use 513:1024 --> (W+1):2*W
            % when k=3 (xz) use 1025:1536 --> (2*W+1):3*W
            % and generally: (k-1)*W+1 : k*W
            ind1 = (k-1)*W + 1 ;
            ind2 = k*W ;
            vid.frames(c).cdata(:,ind1:ind2,1) = img ;
            vid.frames(c).cdata(:,ind1:ind2,2) = img ;
            vid.frames(c).cdata(:,ind1:ind2,3) = img ;
            waitbar(waitbarCounter/Nprogress) ;
        end
    end
    
    % ----------------------------------------------
    %% add time stamp on bottom left corner
    c=0 ;
    prevRealtimeMS = -inf ;
    for it = t:(t+N-1)
        c=c+1 ;
        realtimeMS = round(double(it) / realFrameRate * 1000) ;
        
        if (realtimeMS ~= prevRealtimeMS)
            st = [num2str(realtimeMS) 'ms' ] ;
            if (realtimeMS>0)
                st = ['+' st] ; %#ok<AGROW>
            end
            
            % update timer string
            figure(htimer) ;
            set(htext,'string',st) ;
            
            % save and load image of timer box--should make this better
            print(htimer,'currtime_tempfig','-dpng','-r0') ;
            % saving and loading this image is a time-consuming step.
            % it would be more efficient to
            % have a stack of timer-images pre-made and just use them every
            % time we need a time-stamp.
            
            tm = imread('currtime_tempfig.png') ;
            tm = tm(:,1:70) ;
            %gf = getframe(htimer);
            %tm = gf.cdata(:,1:70) ;
        end
        % embed the timer in the images
        vid.frames(c).cdata(end-Ny+1:end,1:70,1) = tm ;
        vid.frames(c).cdata(end-Ny+1:end,1:70,2) = tm ;
        vid.frames(c).cdata(end-Ny+1:end,1:70,3) = tm ;
        
        % add colored bar for stimulus indicator
        
        if (realtimeMS >=0) && (realtimeMS <= pulseDuration) && ...
                ~isempty(patchColor)
            test_color_bar_im = insertShape(vid.frames(c).cdata,...
                'FilledRectangle', patchCoords ,...
                'Color',{patchColor,patchColor,patchColor});
            vid.frames(c).cdata = test_color_bar_im ;
        end
        
        prevRealtimeMS = realtimeMS ;
    end
    t = it ; % update time
    
    
    % give command line read-out
    disp(['Processing chunk no. ' num2str(chunk) ' / ' num2str(Nchunks)]) ;
    
    % ------------------------------------
    %% write chunk to video file
    for i=1:c
        writeVideo(writerObj,vid.frames(i).cdata);
    end
    
    
    currtime = vid.times(end) + dt ;
    
end
% ----------------------------------------
%% close down video writer and figure
close(h) ;
close(htimer) ;
dos('del currtime_tempfig.png') ;

dur = toc ;
disp (['Overall time for movie ' mp4Filename ' was ' num2str(dur/60) ' mins.']) ;
close (writerObj);
return
end