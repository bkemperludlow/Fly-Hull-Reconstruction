%--------------------------------------------------------------------------
% quick function to show binarized images of fly from 3 views with wings
% and body classified
%--------------------------------------------------------------------------
function h_binIm = showBinarizedImages(analysisOutput, frame_num, h_binIm)

% params and inputs
if ~exist('h_binIm','var')
    h_binIm = figure('position',[ 94   584   560   160]);
end
ww = [-1 1 -1 1] * 48  ;

% read in relevant data
all_fly_bw = analysisOutput.all_fly_bw ; 
body_only_bw = analysisOutput.body_only_bw ; 
CM_pos = analysisOutput.CM_pos ; 

% make plot
set(0, 'CurrentFigure', h_binIm) ;
for cam=1:3
    fly = getImage4D(all_fly_bw,cam,frame_num) ;
    bod = getImage4D(body_only_bw, cam,frame_num) ;
    rgb = zeros(512,512,3,'uint8') ;
    rgb(:,:,1) = uint8(bod)*255 ;
    rgb(:,:,2) = uint8(fly)*255 ;
    subplot(1,3,cam) ;
    imshow(rgb) ;
    xc = squeeze(CM_pos(cam,frame_num,1)) ;
    yc = squeeze(CM_pos(cam,frame_num,2)) ;
    axis([xc xc yc yc]+ww) ;
end
title(frame_num) ;


end