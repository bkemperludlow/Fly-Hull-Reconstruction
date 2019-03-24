function [opened_im,eroded_im] = open_image3D(im,SE)
%Fills the holes of an image, then erodes using the given SE, followed by
%dilation
%Returns both the eroded and opened image.

%Fill up image to correct for uncaptured voxels
im = imfill(im,'holes');
%Then erode away to remove weak connection between wing and elsewhere
eroded_im = imerode(im,SE);
%Identify the components and retain only the largest component
CC=bwconncomp(eroded_im);
[~,wing_idx] = max(cellfun(@numel,CC.PixelIdxList));
eroded_im = zeros(size(eroded_im));
eroded_im(CC.PixelIdxList{wing_idx}) = 1;
%Dilate as the second step of image opening
opened_im = imdilate(eroded_im,SE);
end