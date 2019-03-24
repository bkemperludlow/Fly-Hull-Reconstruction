defineConstantsScript
phantomSDK_setPath ;
LoadPhantomLibraries ;

%dataPath = 'F:\luca\Analysis\pitch up\Expr_7_mov_009\' ;
%datafilename = 'Expr7mov009_Data_manually_corrected.mat' ;
cinPath      = 'F:\luca\07_230514\cines\' ;
%savePath     = 'F:\luca\Analysis\pitch up\Expr_9_mov_009\perturbation movie' ;

exprNumber  = 7 ;
movieNumber = 27 ;

t1 = -86 ; %-38.625 ;
t2 = -15 ; %-15 ;
t3 = 24 ; %13 ;
t4 = 39 ; %32 ;
t5 = 52 ; %55.25 ;
t6 = 60 ; 

imScale = 7 ;
reloadFlag = false ;

startTrackingTime = -900 ;
endTrackingTime = 900 ;
%% Load movies
movieStr    = [ 'Expr ' num2str(exprNumber) ' Movie ' num2str(movieNumber) ' - ' ] ;
%smoothDataFileName = [ outputFilename(1:end-4) '_smooth_data.mat'] ;
%cd(dataPath) 
%if (reloadFlag || ~exist('data','var'))
%    disp(['loading data for expr ' num2str(exprNumber) ' mov ' num2str(movieNumber)]) ;
%    load ([dataPath datafilename]) ;   
%end

%data.params.pulseStartMS = 0 ; 
%data.params.pulseDurationMS = 5.8 ;

%cin_frames = data.params.startTrackingTime : data.params.endTrackingTime ;
cin_frames = startTrackingTime : endTrackingTime ;
fps         = 8000 ;
tsec       =  cin_frames / fps ;
tms        =  tsec * 1000 ;

fr_start = find(tms==t1) ; % will be used to read the images from cin files
fr_end   = find(tms==t6) ;

XZ = 2 ;
XY = 3 ;
YZ = 1 ;

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

%meters_to_pixels = data.params.voxelSize  ;
deg2rad = pi / 180 ; 
Nt = length(tsec) ;

%firsts   = [data.params.metaData.firstImage] ;
%frames   = (0:data.Nimages-1) + data.params.startTrackingTime ;
%frames = (0:length(cin_frames)-1) + startTrackingTime ;
%% Cycle through images
it1 = find(tsec==t1/1000) ;
it2 = find(tsec==t2/1000) ;
it3 = find(tsec==t3/1000) ;
it4 = find(tsec==t4/1000) ;
it5 = find(tsec==t5/1000) ;
it6 = find(tsec==t6/1000) ;

if (isempty(it1) || isempty(it2) || isempty(it3) || isempty(it4) || isempty(it5) || isempty(it6))
    disp('cannot find t in tsec. aborting.') ;
    return ;
end

% load images and project on walls
%imxy = uint8(zeros(512)) ; 
imxzCell = cell(6,1) ; 
imxyCell = cell(6,1) ; 
imyzCell = cell(6,1) ; 
%imyz = uint8(zeros(512)) ;

allBGcell = importdata('allBGcell.mat') ;
%xz_cm = importdata('xz_cm.mat') ;
delta = 38 ;

i = 0 ;
for it = [it1 it2 it3 it4 it5 it6] ;%[it1 it2 it3 it4]
    i = i+1 ;
    currFrame = cin_frames(it) ;
    %{
    cam   = YZ ;
    imyz = uint8(zeros(512)) ;
    imyz_1 = myReadCinImage(cindata{cam}, currFrame) ; % read image
    imyz_2 = imcomplement(imyz_1) - imcomplement(allBGcell{cam,1}) ;   %allBGcell{1,cam}
    
    imyz  = imyz + fliplr(imyz_2) ;
    
    levels_yz = multithresh(imyz,2) ; %3
    temp_yz = imquantize(imyz,levels_yz) ;
    ind1_yz = find(temp_yz == 1) ;
    imyz(ind1_yz) = 0 ;
    imyz_cropped = imyz ; %imcrop(imxz, [512-xz_cm(it, 1)-delta, xz_cm(it, 2)-delta, 2*delta, 2*delta]) ;
    imyzCell{i,1} = imyz_cropped ;
    %}
   
    cam   = XZ ;
    imxz = uint8(zeros(512)) ;
    imxz_1 = myReadCinImage(cindata{cam}, currFrame) ; % read image
    imxz_2 = imcomplement(imxz_1) - imcomplement(allBGcell{1,cam}) ;   %allBGcell{1,cam}
    
    imxz  = imxz + fliplr(imxz_2) ;
    
    levels_xz = multithresh(imxz,2) ;
    temp_xz = imquantize(imxz,levels_xz) ;
    ind1_xz = find(temp_xz == 1) ;
    imxz(ind1_xz) = 0 ;
    imxz_cropped = imxz ; %imcrop(imxz, [512-xz_cm(it, 1)-delta, xz_cm(it, 2)-delta, 2*delta, 2*delta]) ;
    imxzCell{i,1} = imxz_cropped ;
    
    %{
    cam   = data.params.XY ;
    imxy = uint8(zeros(512)) ;
    imxy_1 = myReadCinImage(cindata{cam}, currFrame) ; % read image
    imxy_2 = imcomplement(imxy_1) - imcomplement(allBGcell{1,cam}) ;   %allBGcell{1,cam}
    
    imxy  = imxy + fliplr(imxy_2) ;
    
    levels_xy = multithresh(imxy,3) ;
    temp_xy = imquantize(imxy,levels_xy) ;
    ind1_xy = find(temp_xy == 1) ;
    imxy(ind1_xy) = 0 ;
    imxy_cropped = imxy ; %imcrop(imxz, [512-xz_cm(it, 1)-delta, xz_cm(it, 2)-delta, 2*delta, 2*delta]) ;
    imxyCell{i,1} = imcomplement(imxy_cropped) ;
    %}
end
%htest = fspecial('gaussian') ;

%imxy = imcomplement(imxy) ;
%imxy = imresize(imxy,imScale) ;

%cd('F:\Sam\Pitch Paper\pert diagram')
%cd('F:\Sam\Pitch Paper\MATLAB figs\extreme perturbations')
imSum = uint8(zeros(512)) ;
for j = 2:6
    %figure ; imshow(imxzCell{j,1}) 
    %print(gcf,['xz_' num2str(j)],'-dpng','-r300')
    imSum = imSum + imxzCell{j,1} ;
    %print(gcf,['xy_' num2str(j)],'-dpng','-r300')
end
imSum_scaled = imresize(imSum,imScale) ;
figure ; imshow(imcomplement(imSum_scaled)) ;
%{    
totX = 500 ;
totY = 500 ; 
imxzSum = uint8(zeros(totX,totY)) ;
itVec = [it1 it2 it3 it4 it5] ;

for k = 1:6
    Dist =  int16(2*abs(xz_cm(itVec(k),:) - xz_cm(itVec(1),:))) ;
    if k == 4
        extraDist = [20, 0] ;
    elseif k == 5
        extraDist = [40, 0] ;
    else
        extraDist = [0, 0] ;
    end
    IndX1 = Dist(1) + extraDist(1) + 1 ;
    IndX2 = IndX1+2*delta ;
    IndY2 = totY - Dist(2) + extraDist(2);
    IndY1 = IndY2 - 2*delta ;
    imxzSum(IndY1:IndY2,IndX1:IndX2) = imxzSum(IndY1:IndY2,IndX1:IndX2)...
        +imcomplement(imxzCell{k,1}) ;
    %figure ; imshow(uint8(imxzCell{k,1}))
    %imxz_temp = imresize(uint8(imxzCell{k,1}),imScale) ;
    %print(gcf,['xz_' num2str(k)],'-dpng','-r300')
end
imxzSum = imcomplement(imxzSum) ;
imxzSum_crop = imcrop(imxzSum,[1, totY-225, 400, 225]) ;  
%imyz = imcomplement(imyz) ;
%imyz = imresize(imyz,imScale) ;
%}
%figure ; imshow(imxy)

%figure ; imshow(imyz)
