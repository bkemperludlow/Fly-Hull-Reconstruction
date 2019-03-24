function [wing1_bw, wing2_bw] = findWingsOnly(fly_bw, body_bw, cam)
wing1_bw = [] ;
wing2_bw = [] ;
Nimages = fly_bw.dim(1) ;

for k=1:Nimages
    currfly  = getImage4D(fly_bw, cam, k) ;
    currbody = getImage4D(body_bw, cam, k) ;
    thickbody = bwmorph(currbody,'thicken',1) ;
    dif       = currfly - thickbody ;
    
    figure(88) ;
    subplot(2,3,1) ; imshow(currfly) ;
    subplot(2,3,2) ; imshow(currbody) ;
    subplot(2,3,3) ; imshow(thickbody) ;
    subplot(2,3,4) ; imshow(dif) ;
    pause ;
end