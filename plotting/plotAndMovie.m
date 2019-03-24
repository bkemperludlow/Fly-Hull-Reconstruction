%% Initialize
phantomSDK_setPath ;
LoadPhantomLibraries ;
defineConstantsScript ;

writepath = 'F:\luca\Analysis\pitch up\Expr_7_mov_009\' ;
frameRate = 25 ; 

ExprNum = 7 ;
MovNum = 9 ; 

patchColor = [1 1 1 ] * 0.8 ; 
faceAlpha = 1 ;

%% Get cine file
if (MovNum<10)
    zstr = '00' ;
elseif (MovNum<100)
    zstr = '0' ;
else
    zstr='';
end
cinPath      = 'F:\luca\07_230514\cines\' ;
xzfile = [cinPath 'xz_' zstr num2str(MovNum) '.cin'] ;

cindata = myOpenCinFile(xzfile) ; 

%% Get data

if MovNum < 10
    zstr = '00' ;
elseif MovNum < 100
    zstr = '0' ;
else
    zstr = '' ;
end
datapath = strcat('F:\luca\Analysis\pitch up\Expr_',...
            num2str(ExprNum), '_mov_',zstr,num2str(MovNum)) ;
cd(datapath)
datafilename = strcat(datapath,'\Expr',num2str(ExprNum), ...
        'mov',zstr,num2str(MovNum), '_Data_manually_corrected.mat') ;
load(datafilename) ;        
bodyPitch = data.anglesLabFrame(:,BETA) ;
phiR = -data.anglesBodyFrame(:,PHIR) ; 
phiL = data.anglesBodyFrame(:,PHIL) ; 

%% Create timing things
tstart = -10 ; % in ms
tend   = 30 ; % in ms

pulseStart = 0 ; 
pulseEnd = data.params.pulseLengthMS ;

cin_frames = data.params.startTrackingTime : data.params.endTrackingTime ;
tsec       =  cin_frames / data.params.fps ;
tms        =  tsec * 1000 ;

itstart = find(tsec==tstart/1000) ;
itend   = find(tsec==tend/1000) ;

[~,pulseStartInd] = min(abs(tsec - pulseStart/1000)) ; 
[~,pulseEndInd] = min(abs(tsec - pulseEnd/1000)) ; 

%% Main loop

imxlim = [115 275] ;
imylim = [225 325] ; 
cd([writepath 'flymovie\']) ;

for it=itstart:1:itend
    itstr = num2str(it) ; 
    
    currFrame = cin_frames(it) ;
    imxz = myReadCinImage(cindata, currFrame) ;
    imxzflipped = fliplr(imxz) ; 
    imxzSave = imxzflipped(imylim(1):imylim(2),imxlim(1):imxlim(2)) ;
    %{
    h1 = figure; %('position',[140 50 500 500]) ; 
        set(gcf,'paperPositionMode','auto') ;
        %imshow(imxzflipped,'Border','tight');
        imshow(imxzflipped,'Border','tight');
        set(gca,'ylim',imylim) ;
        set(gca,'xlim',imxlim) ;
    %}    
        
    tempImName = ['flymovie' itstr '.png'] ;
    imwrite(imxzSave, tempImName) ;
    %print(gcf,tempImName,'-dpng','-r300')
    %delete(h1) ; 
end

writerObj1 = VideoWriter('fly_movie.avi') ;
writerObj1.FrameRate = frameRate;
open(writerObj1) ;

for it=itstart:1:itend
    itstr = num2str(it) ;
    tempIm = imread(['flymovie' itstr '.png']) ; 
    writeVideo(writerObj1,tempIm);
end
close(writerObj1)

cd([writepath 'kinematics\']) ;

ylim = [10 180] ; 
xlim = [-10 30] ;
avec = [ylim(1) ylim(1) ylim(2) ylim(2) ylim(1) ] ;

for it=itstart:1:itend
    itstr = num2str(it) ; 
    
    h2 = figure('position',[140 50 1100 500]) ; 
    set(gcf,'paperPositionMode','auto')
    hold on
    set(gca,'fontsize',14)
    set(gca,'ylim',ylim)
    set(gca,'xlim',xlim)
    if it >= pulseStartInd && it<= pulseEndInd
        tsfvec = [tms(pulseStartInd) tms(it) tms(it) tms(pulseStartInd) tms(pulseStartInd) ] ;
        hf = fill(tsfvec , avec,'y') ;
        set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
    elseif it > pulseEndInd
        tsfvec = [tms(pulseStartInd) tms(pulseEndInd) tms(pulseEndInd) tms(pulseStartInd) tms(pulseStartInd) ] ;
        hf = fill(tsfvec , avec,'y') ;
        set(hf,'facecolor',[255 255 153]/255,'facealpha',faceAlpha) ;
    end  
    bp = plot(tms(itstart:it), bodyPitch(itstart:it), 'k.') ;
    rws = plot(tms(itstart:it), phiR(itstart:it), '-ro','LineWidth',1.5,...
        'MarkerSize',5,'MarkerFaceColor','r') ;
    lws = plot(tms(itstart:it), phiL(itstart:it), '-bo','LineWidth',1.5,...
        'MarkerSize',5,'MarkerFaceColor','b') ;
    legend([bp rws lws],{'\theta_B','\phi_R','\phi_L'},'location','northwest')
    xlabel('Time [ms]')
    ylabel('Angle [deg]')
    
    tempImName = ['kinematics' itstr] ; 
    print(gcf,tempImName,'-dpng','-r300')
    delete(h2) ;
    
end

writerObj2 = VideoWriter('kinematics.avi') ;
writerObj2.FrameRate = frameRate ;
open(writerObj2) ;

for it=itstart:1:itend
    itstr = num2str(it) ;
    tempIm2 = imread(['kinematics' itstr '.png']) ; 
    writeVideo(writerObj2,tempIm2);
end
close(writerObj2)


