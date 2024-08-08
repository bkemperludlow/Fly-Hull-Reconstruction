% -------------------------------------------------------------------------
% create triptych video from three cine files. this version ("*Sam.m")
% allows for unequal video sizes
% -------------------------------------------------------------------------
function stat = combineMoviesFromCinSam(metaData, outputFileName,fps, ...
    realFrameRate)

tic ;

stat = 0 ;

XZ = 1 ;
XY = 2 ;
YZ = 3 ;

dt = 1/fps ;

chunksize = 365 ;

if (~exist('realFrameRate','var'))
    realFrameRate = 8000 ;
end

mp4Filename = [outputFileName '.mp4'] ;
conf.videoQuality = 75;

% stuff for embedding the timer on the bottom left corner
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
htext = text( 'units','pixels','position',[1 0],'fontsize',12,...
    'FontWeight','bold', 'VerticalAlignment','bottom','string',st) ;

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

movObj = cell(1,3) ;

cams = [XZ XY YZ] ;
availableCams = cams(flags) ;

Hvec = [metaData.height] ;
Wvec = [metaData.width] ;

Hvec = unique(Hvec(flags)) ;
Wvec = unique(Wvec(flags)) ;

H = max(Hvec) ; W = max(Wvec) ; clear Hvec Wvec ;

Nchunks = ceil(Ntot/chunksize) ;

t = globalFirstImage ; % current "time" counted in frames
currtime = 0 ;

ttl = ['Processing movie ' mp4Filename ] ;
ttl(ttl=='_')=' ' ;

h = waitbar(0,ttl);

waitbarCounter = 0 ;

writerObj=VideoWriter(mp4Filename,'MPEG-4');

open(writerObj);


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
                %img = read(obj,n) ;
                %[img, ~] = ReadCineFileImage(metaData(k).filename, it,false); % old code with memory leak
                %% modified for now to get brighter images. temporary fix 
                img = myReadCinImage(metaData(k).cindata, it) +20;
            catch ERR
                disp(ERR) ;
                keyboard;
            end
            waitbarCounter = waitbarCounter + 1 ;
            % adjsut image dimensions to match largest frame
            [sz1, sz2] = size(img) ;
            if (sz1 ~= H)
                pad_height = round((H - sz1)/2) ;
                img = padarray(img, [pad_height, 0], 0, 'both') ;
            end
            if (sz2 ~= W)
                pad_width = round((W - sz2)/2) ;
                img = padarray(img, [0, pad_width], 0, 'both') ;
            end
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
    
    %     keyboard ;
    
    % add time stamp on bottom left corner
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
            
            figure(htimer) ;
            
            % old code
            % imshow(I);
            % set(gca,'position',[0 0 1 1 ]) ;
            % set(gca,'units','pixels') ;
            % set(gca,'position',[0 0 Nx Ny]) ;
            % set(gcf,'position',[500 700 Nx Ny]) ;
            % set(gca,'position',[0 0 Nx Ny]) ;
            set(htext,'string',st) ;
            %htext = text( 'units','pixels','position',[1 0],'fontsize',12,'FontWeight','bold',...
            %    'VerticalAlignment','bottom','string',st) ;
            
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
        prevRealtimeMS = realtimeMS ;
    end
    t = it ; % update time
    
    
    % add vid to movie
    
    disp(['Processing chunk no. ' num2str(chunk) ' / ' num2str(Nchunks)]) ;
    
    
    %
    %     if (chunk==1)
    %         %mmwrite(wmvFilename,vid,conf,'Continue')
    %
    %     elseif (chunk<Nchunks) % middle chunk
    %         %mmwrite(wmvFilename,vid,conf,'Continue','Initialized')
    %
    %     else % last chunk
    %         %mmwrite(wmvFilename,vid,conf,'Initialized')
    %
    %     end
    %
    %     if (chunk<Nchunks)
    %         %disp('Pausing 5 sec') ;
    %         %pause(5);
    %     end
    
    for i=1:c
        writeVideo(writerObj,vid.frames(i).cdata);
    end
    
    
    currtime = vid.times(end) + dt ;
    
end

close(h) ;
close(htimer) ;
dos('del currtime_tempfig.png') ;

dur = toc ;
disp (['Overall time for movie ' mp4Filename ' was ' num2str(dur/60) ' mins.']) ;
close (writerObj);
return