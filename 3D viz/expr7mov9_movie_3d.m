% for k=1:1203 ; im=imread(['Expr22Mov4_v2_' num2str(k) '.png']) ; imwrite(im,'test.tif','writemode','append') ; disp(k) ; end ;
%{
metaData.firstImage = 1 ;
metaData.lastImage  = 1203 ;
metaData.filename   = 'test' ;
metaData.height     = 768 ;
metaData.width      = 1024 ;

%}
defineConstantsScript
phantomSDK_setPath ;
LoadPhantomLibraries ;

tstart = -15 ; % in ms
tend   = 40 ; % +43.5; % in ms

saveImagesFlag =  0 ;

datafilename = 'F:\luca\Analysis\pitch down\Expr_7_mov_065\Expr7mov065_Data_manually_corrected.mat' ;
cinPath      = 'F:\luca\07_230514\cines\' ;
savePath     = 'F:\luca\Analysis\pitch down\Expr_7_mov_065\perturbation movie' ;

exprNumber  = 7 ;
movieNumber = 65 ;
Lpix = 512 ;
imageResolution = [1024 768] ; % [800 600] ;

%axlim = [ 50 315 150 350 125 315 ] ; for expr7mov9
axlim = [ 50 315 100 350 125 315 ] ; % for expr7mov65
%axlim = [ 30 (Lpix-100) 1 Lpix 180 (Lpix-40) ] ; % for expr22mov4 entire movie

%axlim = [ 100-10 (Lpix-180) 170 (Lpix-20) 210 (Lpix-100) ] ; % for expr22mov4 fig 1a (smaller volume)
scale = 10 ;
%axlim = [ 40 Lpix 100 Lpix-20 200 350] ;

az = 40 ; %54 % used in fig1a
el = 14 ; % 12 ; 


% ----------------------------------------------------
% BASED ON The SCRIPT flyPlotMain1a.m
% AND THE CROPPING DONE IN THE SCRIPT figure_1a_v2.m
% 
% see a similar movie without rotating the camera in
% movie1v3_script.m
% ----------------------------------------------------

% change lighting such that bg is gray
% superpose images and 3D-flies
% maybe add sticks from the floor to the flies/traj
% fly's material and light source

% zoom and cam distance
%{
camera {
    location <-300, -100, 500> 
    angle 55         // FOV angle
    right -4/3*x     // don't touch - defines a right handed coordinate system
    sky <0,0,1>      // which way is up
    look_at <128, 159, 175>   // center of FOV
}

%}
% campos camtarget camlookat

%  set(gca,'cameraViewAngle',55) ; campos([-300 -100 350]) ; camtarget([200 256 100])

format compact

reloadFlag = false ;
%{
PC = getPC() ; % lab PC
switch (PC)
    case 1
        tifPath = 'F:\Roll\Data\22_090211\' ;
    case 2
        tifPath = NaN ; % 'C:\Tsevi\Roll\Data\22_090211\' ;
end
%}

%% LOAD DATA

%outputFilename = getDataFileName(exprNumber, movieNumber, PC) ;
movieStr    = [ 'Expr ' num2str(exprNumber) ' Movie ' num2str(movieNumber) ' - ' ] ;
%smoothDataFileName = [ outputFilename(1:end-4) '_smooth_data.mat'] ;

if (reloadFlag || ~exist('data','var'))
    disp(['loading data for expr ' num2str(exprNumber) ' mov ' num2str(movieNumber)]) ;
    load (datafilename) ;
end

data.params.pulseStartMS = 0 ; 
data.params.pulseDurationMS = 5.8 ;

cin_frames = data.params.startTrackingTime : data.params.endTrackingTime ;
tsec       =  cin_frames / data.params.fps ;
tms        =  tsec * 1000 ;

fr_start = find(tms==tstart) ; % will be used to read the images from cin files
fr_end   = find(tms==tend) ;

XZ = data.params.XZ ;
XY = data.params.XY ;
YZ = data.params.YZ ;

% open cin files

if (movieNumber<10)
    zstr = '00' ;
elseif (movieNumber<100)
    zstr = '0' ;
else
    zstr='';
end

xyfile = [cinPath 'xy_' zstr num2str(movieNumber) '.cin'] ;
xzfile = [cinPath 'xz_' zstr num2str(movieNumber) '.cin'] ;
yzfile = [cinPath 'yz_' zstr num2str(movieNumber) '.cin'] ;

cindata = cell(1,3) ;
cindata{XY} = myOpenCinFile(xyfile) ; 
cindata{YZ} = myOpenCinFile(yzfile) ; 
cindata{XZ} = myOpenCinFile(xzfile) ; 

% ------------------------------

meters_to_pixels = data.params.voxelSize  ;
deg2rad = pi / 180 ; 
Nt = length(tsec) ;

%firsts   = [data.params.metaData.firstImage] ;
frames   = (0:data.Nimages-1) + data.params.startTrackingTime ;

% figure ;
% h1 = surface(x,y,z) ; % z=0 plane
% h2 = surface(z,x,y) ; % x=0 plane
% h3 = surface(x,z,y) ; % y=0 plane

Lpix = 512 ;
%Lpix = size(imxy,1) ;
%[x,y] = meshgrid([1 Lpix]) ;
%z = zeros(size(x)) ;

%axlim = [ 80 Lpix 1 256 128 328 ] ;


[x1,y1] = meshgrid(axlim(1:2), axlim(3:4));
z1 = zeros(size(x1)) + axlim(5);

[y2, z2] = meshgrid(axlim(3:4), axlim(5:6));
x2 = zeros(size(y2)) + axlim(2); 

y2b = y2 ;
z2b = z2 ;
x2b = zeros(size(y2)) + axlim(1); 

[x3, z3] = meshgrid(axlim(1:2), axlim(5:6));
y3 = zeros(size(x3)) + axlim(3); 

x3b = x3 ;
z3b = z3 ;
y3b = zeros(size(x3)) + axlim(4); 

% evaluate smoothed coordinates
% -----------------------------

% use tsec as the time vector, in seconds (not ms)
% later need to find the frame corresponding to a give time
% and to plot perturbation (with the correct time offset)
% use real units for the axis

defineConstantsScript ;

%rollEstErr = 1;
%[anglesLabFrame, anglesBodyFrame, t, newEtaLab, newEtaBody,sp_rho, smoothed_rho, rho_t, rho_samp] ... 
%    = calcAngles_quick_and_dirty(data, rollEstErr, true) ;

xy_cm = importdata('xy_cm.mat') ;
xz_cm = importdata('xz_cm.mat') ;
yz_cm = importdata('yz_cm.mat') ;
x_cm_avg = mean([512-xy_cm(:,1), 512-xz_cm(:,1)], 2) ;
y_cm_avg = mean([xy_cm(:,2), 512-yz_cm(:,1)], 2) ;
z_cm_avg = mean([512-xz_cm(:,2), 512-yz_cm(:,2)], 2) ;

bodyRcm  = [x_cm_avg, y_cm_avg, z_cm_avg] ;
%CMEstErr = 3 ; 
%[~,bodyRcm,~] = mySplineSmooth(tsec,bodyRcm_raw, CMEstErr);
%bodyRcm = bodyRcm' ; 

bodyyaw_raw = data.anglesLabFrame(:, PHIB) ;
bodypitch_raw  = data.anglesLabFrame(:, THETAB) ;
bodyroll_raw  = data.anglesLabFrame(:, PSIB) ;

wingstrokeR_raw = data.anglesBodyFrame(:,PHIR) ;
wingdevR_raw =  data.anglesBodyFrame(:,THETAR) ; 
wingpitchR_raw = data.anglesBodyFrame(:,ETAR) ; 

wingstrokeL_raw = data.anglesBodyFrame(:,PHIL) ; 
wingdevL_raw =  data.anglesBodyFrame(:,THETAL) ; 
wingpitchL_raw = data.anglesBodyFrame(:,ETAL) ; 

%bodyYPR_raw  = data.anglesLabFrame(:, [PHIB THETAB PSIB]) ;
%rightYPR_raw  = data.anglesLabFrame(:, [PHIR THETAR ETAR]) ;
%leftYPR_raw  = data.anglesLabFrame(:, [PHIL THETAL ETAL]) ;

bodyyawEstErr = .5 ; %deg
bodypitchEstErr = .5 ; 
bodyrollEstErr = .5 ; 
wingstrokeEstErr = 2 ; 
wingdevEstErr = 3 ; 
wingpitchEstErr = 5 ; 

[~, bodyyaw, ~] = mySplineSmooth(tsec,bodyyaw_raw, bodyyawEstErr);
[~, bodypitch, ~] = mySplineSmooth(tsec,bodypitch_raw, bodypitchEstErr);
[~, bodyroll, ~] = mySplineSmooth(tsec,bodyroll_raw, bodyrollEstErr);

[~, wingstrokeR, ~] = mySplineSmooth(tsec,wingstrokeR_raw, wingstrokeEstErr);
[~, wingdevR, ~] = mySplineSmooth(tsec,wingdevR_raw, wingdevEstErr);
[~, wingpitchR, ~] = mySplineSmooth(tsec,wingpitchR_raw, wingpitchEstErr);

[~, wingstrokeL, ~] = mySplineSmooth(tsec,wingstrokeL_raw, wingstrokeEstErr);
[~, wingdevL, ~] = mySplineSmooth(tsec,wingdevL_raw, wingdevEstErr);
[~, wingpitchL, ~] = mySplineSmooth(tsec,wingpitchL_raw, wingpitchEstErr);


bodyYPR = [bodyyaw' bodypitch' bodyroll'] ;
rightYPR = [wingstrokeR' wingdevR' wingpitchR'] ;
leftYPR = [wingstrokeL' wingdevL' wingpitchL'] ;

% unit conversions
%bodyRcm  = bodyRcm  + data.params.volLength / (data.params.voxelSize) ; % * meters_to_pixels  + data.params.volLength/meters_to_pixels/2 ;
bodyYPR  = bodyYPR  * deg2rad ;
rightYPR = rightYPR * deg2rad ;
rightYPR(:,1) = -rightYPR(:,1) ;
leftYPR  = leftYPR  * deg2rad ;

resolution =  20 ;

%%minr = min(bodyRcm) ;
%%maxr = max(bodyRcm) ;
%axlim = [minr(1) maxr(1) minr(2) maxr(2) minr(3) maxr(3) ] + ...
%    [-1 1 -1 1 -1 1]*40 ;
%axlim = [ 1 Lpix 1 Lpix 1 Lpix ] ;

%%clear minr maxr ;

mainfig = figure('position',[500  100   imageResolution],'units','pixels') ; % size(imageResolution)=[1,2] ;
ax = gca ; % will be deleted and recreated later to make htimer on front
axpos = get(ax,'position') ;

axtimer = axes('units','pixels') ;
delete(ax) ;
ax = axes('parent',mainfig,'position',axpos) ;
clear axpos ;

set(axtimer,'position',[ 1 1 45 12],'xtick',[], 'ytick',[],'xcolor',[1 1 1],'ycolor',[1 1 1]) ;
timerTextPos=[0 0.5] ;

axes(axtimer)
htext = text(timerTextPos(1),timerTextPos(2),' Timer','fontweight','bold') ;

set(gcf,'paperpositionMode','auto') ;
set(gcf,'color',[1 1 1 ]);% [16 16 16 ] / 255) ;
set(gcf,'invertHardcopy','off')

set(gcf,'renderer','zbuffer')



pinType = 2 ; % 0==no pin, 1==roll pin, 2==pitch
thetab0 = 45 * pi / 180 ;

[flyGrp, bodyGrp, rightWingGrp, leftWingGrp, dL, rightWingTip, leftWingTip] ...
    = draw3Dfly_v3(ax, scale, resolution, pinType, thetab0); %#ok<ASGLU>

box off ; grid on ; view(3) ; rotate3d on ; axis equal ;

axis auto ;
dgrid = 50 ;

xtick = [];%axlim(1):dgrid:axlim(2) ;
ytick = [];%axlim(3):dgrid:axlim(4) ;
ztick = [];%axlim(5):dgrid:axlim(6) ;

set(gca,'xtick',xtick,'ytick',ytick,'ztick',ztick) ;
set(gca,'xticklabel',[],'yticklabel',[],'zticklabel',[]) ;
set(gca,'Ticklength' , [0 0]) ;

camproj('orthographic') ;%camproj('perspective')

pulseStartTime = data.params.pulseStartMS / 1000 ;
pulseEndTime   = pulseStartTime + data.params.pulseDurationMS / 1000 ;

pulseStartInd = find(tsec==pulseStartTime) ;
tempDiff = abs(tsec - pulseEndTime) ;  
[~,pulseEndInd]   = min(tempDiff) ;

pulseInd = pulseStartInd:pulseEndInd ;
hcm1 = [] ;
hcm2 = [] ;
hcm3 = [] ;
hcm4 = [] ;
hrcontrail = [] ;
hlcontrail = [] ;

rightTipPos = zeros(Nt,3) ;
leftTipPos  = zeros(Nt,3) ;

imgCounter = 0 ;

% "walls" for the three images
h1  = surface(x1,y1,z1) ; % z=0 plane
h2  = surface(x2,y2,z2) ; % x=0 plane
h2b = surface(x2b,y2b,z2b) ; % x=0 plane
h3  = surface(x3,y3,z3) ; % y=0 plane % old: surface(x,z+Lpix,y)
h3b = surface(x3b,y3b,z3b) ; 

% set lighting properties
hlight = light ;
set(hlight,'position',[-0.928136 -1.27747 0.711783]) ;

material shiny ;
set([h1 h2 h3 h2b h3b], 'ambientStrength',0.4);  
set([h2b h3b], 'ambientStrength',0.9);  

% set lighting properties for all surfaces in the fly
c = findobj(flyGrp,'Type','surface');
set(c,'ambientstrength',0.7); % note that this line comes after "material shiny"



cycle = 0 ; 
contrailDT = 18 ;

% find initial and final indiced that correspond to tstart and tend
itstart = find(tsec==tstart/1000) ;
itend   = find(tsec==tend/1000) ;
if (isempty(itstart))
    disp('cannot find tstart in tsec. aborting.') ;
    return ;
end

if (isempty(itend))
    disp('cannot find tend in tsec. aborting.') ;
    return ;
end

if saveImagesFlag
    cd(savePath) ;
end

while(cycle<=2) 
    cycle=cycle+1 ;
  
    %for it=1:1:length(tsec)
    %currFrameInd = fr_start - 1 ;
    for it=itstart:1:itend
        %currFrameInd = currFrameInd + 1 ;
        
        delete(hcm1) ;
        delete(hcm2) ;
        try
            delete(hcm3) ;
        catch
            % do nothing ;
        end
        try
            delete(hcm4) ;
        catch
            % do nothing ;
        end
        delete(hrcontrail) ;
        delete(hlcontrail) ;
        
        setFlyDOF(flyGrp, rightWingGrp, leftWingGrp, ...
            bodyRcm(it,:), bodyYPR(it,:), rightYPR(it,:) , leftYPR(it,:), thetab0, dL )
        
        MB = get(flyGrp,'Matrix') ;
        MR = get(rightWingGrp,'Matrix') ;
        ML = get(leftWingGrp,'Matrix') ;
        
        rr = [rightWingTip ; 1] ;
        r = (MB * MR * rr)' ;
        rightTipPos(it,:) = r(1:3) ;
        
        rl = [leftWingTip ; 1] ;
        r = (MB * ML * rl)' ;
        leftTipPos(it,:) = r(1:3) ;
        
        axis equal ; 
        axis(axlim);
        %axis([1 512 1 512 1 512]);
        
        hold on ;
        if (it<pulseStartInd)
            greenInd = itstart:it ;
        elseif (it>=pulseStartInd && it<=pulseEndInd)
            greenInd = itstart:(pulseStartInd-1) ;
        elseif (it>pulseEndInd)
            greenInd = [itstart:(pulseStartInd-1) (pulseEndInd+1):it] ;
        end
        
        temprcm = NaN(size(bodyRcm)) ;
        temprcm(greenInd,:) = bodyRcm(greenInd,:) ;
        
        %hcm1 = plot3(bodyRcm(greenInd,1), bodyRcm(greenInd,2),bodyRcm(greenInd,3),'g-','linewidth',2) ;
        %hcm2 = plot3(bodyRcm(greenInd,1), bodyRcm(greenInd,2),zeros(1,length(greenInd))+axlim(5),'-','color',[1 1 1],'linewidth',2) ;
        tmpind = itstart:it ;
        hcm1 = plot3(temprcm(tmpind,1), temprcm(tmpind,2),temprcm(tmpind,3),'-','linewidth',2,'color',[0 0.75 0.1]) ;
        hcm2 = plot3(temprcm(tmpind,1), temprcm(tmpind,2),zeros(1, length(tmpind))+axlim(5),'-','color',[1 1 1],'linewidth',2) ;
        clear tmpind ;
        
        if (it>=pulseStartInd)
            it2 = min([it pulseEndInd]) ;
            hcm3 = plot3(bodyRcm(pulseStartInd:it2,1), bodyRcm(pulseStartInd:it2,2),bodyRcm(pulseStartInd:it2,3),...
                'r-','linewidth',3) ;
            nz = it2-pulseStartInd+1 ;
            hcm4 = plot3(bodyRcm(pulseStartInd:it2,1), bodyRcm(pulseStartInd:it2,2),zeros(1,nz)+axlim(5),...
                '-','color',[0.6 0.05 0.05],'linewidth',3) ;
        end
        
        if (cycle>2) ;
            tc = max([itstart (it-contrailDT-(cycle>3)*contrailDT*4)]) ;
            
            hrcontrail = plot3(rightTipPos(tc:it,1), rightTipPos(tc:it,2), rightTipPos(tc:it,3), ...
                '.-','color',[0.6 0.0 0.0]) ;
            
            hlcontrail = plot3(leftTipPos(tc:it,1), leftTipPos(tc:it,2), leftTipPos(tc:it,3), ...
                'b.-') ;
        end
        %title(it) ;
        %plot3(leftTipPos(it,1),leftTipPos(it,2),leftTipPos(it,3),'bs','markersize',15) ;
        %view(az, 14) ;
        view(az,el) ;
        if (cycle>1 && az<111)
            az = az + 0.1 ;
        end
        if (cycle>1 && az<111)
            d = 3/1000 ;
            el = el - d ;
        end
            
        % load images and project on walls
        cam   = data.params.YZ ; %n = frames(it) - firsts(cam) + 1 ;
        currFrame = cin_frames(it) ;
        
        imyz = myReadCinImage(cindata{cam}, currFrame) ; % read image
        %imyz  = imread(imageFileNames{cam},'index', n,'info',infoCell{cam}) ; % read image
        imyz  = fliplr(imyz) ;
        
        cam   = data.params.XZ ; %n = frames(it) - firsts(cam) + 1 ;
        imxz = myReadCinImage(cindata{cam}, currFrame) ; % read image
        %imxz  = imread(imageFileNames{cam},'index', n,'info',infoCell{cam}) ; % read image
        imxz  = fliplr(imxz) ;
        
        cam   = data.params.XY ; %n = frames(it) - firsts(cam) + 1 ;
        imxy = myReadCinImage(cindata{cam}, currFrame) ; % read image
        %imxy  = imread(imageFileNames{cam},'index', n,'info',infoCell{cam}) ; % read image
        imxy  = rot90(imxy,2) ;
        
        % crop images - code taken from figure_1a_v2.m
        % --------------------------------------------
        im1 = flipud(double(imxy)) ;
        im1 = im1(axlim(3):axlim(4), axlim(1):axlim(2)) ;
        
        im2 = flipud(double(imyz)) ;
        im2 = im2(axlim(5):axlim(6), axlim(3):axlim(4)) ;
        
        im3 = flipud(double(imxz)) ;
        im3 = im3(axlim(5):axlim(6), axlim(1):axlim(2)) ;
        % --------------------------------------------
        
        set(h1,'CData',(double(im1)),'FaceColor','texturemap','linestyle','-')
        set(h2,'CData',(double(im2)),'FaceColor','texturemap','linestyle','-')
        set(h2b,'CData',(double(im2)),'FaceColor','texturemap','linestyle','-')
        set(h3,'CData',(double(im3)),'FaceColor','texturemap','linestyle','-')
        set(h3b,'CData',(double(im3)),'FaceColor','texturemap','linestyle','-')
        
        colormap gray ;
        
        if (az<0)
            set(h2,'visible','on') ;
            set(h2b,'visible','off') ;
        else
            set(h2,'visible','off') ;
            set(h2b,'visible','on') ;
        end
        
        if (az<=-90) || (az>=90)
            set(h3,'visible','on') ;
            set(h3b,'visible','off') ;
        else
            set(h3,'visible','off') ;
            set(h3b,'visible','on') ;
        end
        
        set(gca,'position',[0.01 0.01 0.98 0.99]) ;
        
        % code from flyPlotMain1a.m
        %set(gca,'cameraViewAngle',5.7270 );
        %camtarget([244.4331  222.5500  182.0753]) ;
        %campos([-3.5505   -2.6371    2.1019]*1000) ;
        
        % EMBED TIMER in ms
        timeVal = floor(tsec(it)*1000) ;
        if (timeVal<0)
            timeStr = [ ' ' num2str(timeVal) 'ms' ];
        else
            timeStr = [ '  ' num2str(timeVal) 'ms' ];
        end
        
        set(htext,'String',timeStr) ;
        
        if (it>=pulseStartInd && it<=pulseEndInd)
            set(htext,'color','r') ;
        else
            set(htext,'color','k') ;
        end
        
        view(az,el) ;
        imgCounter = imgCounter + 1;
        if (saveImagesFlag)
            print(gcf,['perturbation_animation_' num2str(imgCounter)],'-dpng','-r300') ;
            % '-r0' below prints in the screen resolution imageResolution (see above)
            %print(gcf,['.\images\Expr' num2str(exprNumber) 'Mov' num2str(movieNumber) '_' num2str(imgCounter)],'-dpng','-r0') ;
        end
        disp(imgCounter)
        pause(0.01) ;
        %keyboard ;
        %return

    end
end
clear nz timeStr timeVal

% close cin files
for cam=1:3
    myCloseCinFile(cindata{cam}) ;
end

%{
while(1)
    for t=linspace(0,2*pi,100)
        phi = 90 + 90 * cos(t) ;
        psi = 90 - 45 * sin(t) ;
        rightYPR = [phi 0 psi] * pi/180;
        leftYPR  = [phi 0 psi] * pi/180;
        bodyYPR  = [30 thetab0*180/pi 0 ] * pi / 180;
        bodyRcm  = [3 3 1 ] ;
        
        setFlyDOF(flyGrp, rightWingGrp, leftWingGrp, ...
            bodyRcm, bodyYPR, rightYPR , leftYPR, thetab0, dL )
        axis([-5 10 -5 10 -4 10]) ;
        
        pause(0.01) ;
        
    end
end
%}