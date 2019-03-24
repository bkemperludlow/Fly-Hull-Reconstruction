function [x, y] = inputPointsFromCine(fileName, frameIndices)
% before calling this function call:
% phantomSDK_setPath ;
% LoadPhantomLibraries();
%
% a single image can be read simply by: 
% [currim, ~] = ReadCineFileImage(fileName, frameNumber, false);
% no need to use all the other crap

N = length(frameIndices) ;

% get cine parameters (copied from ReadCineFileImage)
% original function is in: F:\Tsevi\MatlabCode\PhantomSDK\PhMatlabSDK\SimpleCineFileReader
% check which copy of the package needs to be used.

% LoadPhantomLibraries() is in: % F:\Tsevi\MatlabCode\PhantomSDK\PhMatlabSDK\PhMatlab\Other

%% Create the cine handle from the cine file.
%Is recomended that cine handle creation should be done once for a batch of image readings. 
%This will increase speed.
[HRES, cineHandle] = PhNewCineFromFile(fileName);
if (HRES<0)
	[message] = PhGetErrorMessage( HRES );
    error(['Cine handle creation error: ' message]);
end
%% read the saved range
pFirstIm = libpointer('int32Ptr',0);
PhGetCineInfo(cineHandle, PhFileConst.GCI_FIRSTIMAGENO, pFirstIm);
firstIm = pFirstIm.Value;
pImCount = libpointer('uint32Ptr',0);
PhGetCineInfo(cineHandle, PhFileConst.GCI_IMAGECOUNT, pImCount);
lastIm = int32(double(firstIm) + double(pImCount.Value) - 1);

if (min(frameIndices)<firstIm || max(min(frameIndices))>lastIm)
    error(['Image number must be in saved cine range [' num2str(firstIm) ';' num2str(lastIm) ']']);
end

% get image width and height - see "Phantom SDK Reference Manual.pdf" page 99
pWidth = libpointer('int32Ptr',0);
PhGetCineInfo(cineHandle, PhFileConst.GCI_IMWIDTH, pWidth);
width = pWidth.Value;

pHeight = libpointer('int32Ptr',0);
PhGetCineInfo(cineHandle, PhFileConst.GCI_IMHEIGHT, pHeight);
height = pHeight.Value;

%% pre-allocate memory (need to change for 16bit if needed)
ims = zeros(N, height, width,'uint8') ;
x   = zeros(N,1) ;
y   = zeros(N,1) ;

%% read all images
for k=1:N
    [currim, ~] = ReadCineFileImage(fileName, frameIndices(k), false);
    ims(k,:,:) = currim ;
end
clear currim k

winx = 50 ;
winy = 50 ;

% get first point
h1=figure ;
imshow(squeeze(ims(1,:,:))) ;
title('Click generally where the fly center is') ;
[x0, y0] = ginput(1) ;
axvec   = [x0-winx x0+winx y0-winy y0+winy] ;
close(h1) ;

h2 = figure ;
for k=1:N
    imshow(squeeze(ims(k,:,:))) ;
    title(frameIndices(k)) ;
    axis(axvec) ;
    [x0, y0] = ginput(1) ;
    x(k) = x0 ; 
    y(k) = y0 ;
    axvec   = [x0-winx x0+winx y0-winy y0+winy] ;
    pause(0.1) ;
end

close (h2) ;

figure ; plot(x,y,'ko-','markerfacecolor','y') ; 
axis equal ; box on ; grid on ;

disp('done.')
end
