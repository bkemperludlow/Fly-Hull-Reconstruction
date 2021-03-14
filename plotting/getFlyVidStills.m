%--------------------------------------------------------------------------
% function to grab stills from a given movie for all three cameras
%
% INPUTS:
%   -cinRoot: path to folder containing movie .cin files (e.g. 'D:\Fly
%       Data\VNC MN Chrimson\04_10112017\') [string]
%   -ExprNum: experiment number 
%   -MovNum: movie number to take images from
%   -snapshotTimes: array of time values (in seconds) for the desired
%       stills
%   -scale: scalar value to resize image by.
%   -trimFlag: boolean value. if true, trim images to small region around
%       the fly. if false, take full image (typically 512x512)
%   -plotFlag: boolean value. Show image output?
%   -saveFlag: boolean value. Save results?
%   -savePath: path to folder where results should be saved
%
% OUTPUTS:
%   -im_struct_*: data structure containing output images. * = camera name 
%       (xy, xz, yz)
%   -im_*_final: output image with all snapshots merged.
%       
%
% To further process these images, I'll typically:
%   -Open all images as layers in GIMP
%   -Turn white channel to alpha (there exists an add-on to do this across
%   layers). NB: Need to select Image->Mode->RGB to alter alpha
%   -Create white background layer
%   -pick times to show, hide the rest
%   -export to png
%
%--------------------------------------------------------------------------
function [im_struct_xy, im_struct_xz, im_struct_yz, im_xy_final, ...
    im_xz_final, im_yz_final] = getFlyVidStills(cinRoot, ExprNum, MovNum, ...
     snapshotTimes, scale, trimFlag, plotFlag, saveFlag, savePath)
% ---------------------------------
%% params and inputs
if ~exist('scale','var') || isempty(scale)
    scale = 5.0 ; % scale used to resize images
end
if ~exist('trimFlag','var') || isempty(trimFlag)
    trimFlag = false ; %trim images to be mostly fly?
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = false ; %plot combined image?
end
if ~exist('saveFlag','var') || isempty(saveFlag)
    saveFlag = false ; 
end
if ~exist('savePath','var') || isempty(savePath)
    savePath = [] ; 
    saveFlag = false ; 
end

% -------------------------------------------------------------------------
%% get flight stills stored in structures
im_struct_xy = ...
    get_fly_snapshots(cinRoot, ExprNum, MovNum, 'XY',...
     snapshotTimes,savePath, saveFlag, trimFlag, scale)  ;
im_struct_xz = ...
    get_fly_snapshots(cinRoot, ExprNum, MovNum, 'XZ',...
     snapshotTimes,savePath, saveFlag, trimFlag, scale)  ;
im_struct_yz = ...
    get_fly_snapshots(cinRoot, ExprNum, MovNum, 'YZ',...
     snapshotTimes,savePath, saveFlag, trimFlag, scale)  ;

% -------------------------------------------------------------------------
%% save image structures?
if saveFlag
    save(fullfile(savePath, 'im_struct_xy.mat'), 'im_struct_xy')
    save(fullfile(savePath, 'im_struct_xz.mat'), 'im_struct_xz')
    save(fullfile(savePath, 'im_struct_yz.mat'), 'im_struct_yz')
end

% -------------------------------------------------------------------------
%% create overlay images
im_xy_comb = uint8(zeros(size(im_struct_xy(1).image))) ;
im_xz_comb = uint8(zeros(size(im_struct_xz(1).image))) ;
im_yz_comb = uint8(zeros(size(im_struct_yz(1).image))) ;

loop_range = 1:length(snapshotTimes) ;

for i = loop_range
    im_xy = im_struct_xy(i).image ;
    im_xz = im_struct_xz(i).image ;
    im_yz = im_struct_yz(i).image ;
    
    % ---------------------
    % xy image
    min_xy_val = double(min(im_xy(:)))/255.0 ;
    im_xy_adjust = imadjust(im_xy, [min_xy_val, 1.0],[]) ;
    im_xy_complement = imcomplement(im_xy_adjust) ;
    im_xy_comb = imadd(im_xy_comb, im_xy_complement) ;
    im_struct_xy(i).image_processed = im_xy_adjust ;
    
    % ---------------------
    % xz image
    min_xz_val = double(min(im_xz(:)))/255.0 ;
    im_xz_adjust = imadjust(im_xz, [min_xz_val, 1.0],[]) ;
    im_xz_complement = imcomplement(im_xz_adjust) ;
    im_xz_comb = imadd(im_xz_comb, im_xz_complement) ;
    im_struct_xz(i).image_processed = im_xz_adjust ;
    
    % ---------------------
    % yz image
    min_yz_val = double(min(im_yz(:)))/255.0 ;
    im_yz_adjust = imadjust(im_yz, [min_yz_val, 1.0],[]) ;
    im_yz_complement = imcomplement(im_yz_adjust) ;
    im_yz_comb = imadd(im_yz_comb, im_yz_complement) ;
    im_struct_yz(i).image_processed = im_yz_adjust ;
    
    % save processed images
    if saveFlag
        imwrite(im_xy_adjust, fullfile(savePath, ['Expr_' num2str(ExprNum) '_mov_' ...
            num2str(MovNum,'%03d') '_XY_' num2str(im_struct_xy(i).frame,'%04d')...
            '_processed.png'])) ;
        imwrite(im_xz_adjust,  fullfile(savePath, ['Expr_' num2str(ExprNum) '_mov_' ...
            num2str(MovNum,'%03d') '_XZ_' num2str(im_struct_xz(i).frame,'%04d')...
            '_processed.png'])) ;
        imwrite(im_yz_adjust,  fullfile(savePath, ['Expr_' num2str(ExprNum) '_mov_' ...
            num2str(MovNum,'%03d') '_YZ_' num2str(im_struct_yz(i).frame,'%04d')...
            '_processed.png'])) ;
    end
end

% -------------------------------------------------------------------------
%% take complement of image and show results
im_xy_final = imcomplement(im_xy_comb) ;
im_xz_final = imcomplement(im_xz_comb) ;
im_yz_final = imcomplement(im_yz_comb) ;

if plotFlag
    figure ;
    subplot(1,3,1)
    imshow(im_xz_final)
    title('XZ')
    
    subplot(1,3,2)
    imshow(im_xy_final)
    title('XY')
    
    subplot(1,3,3)
    imshow(im_yz_final)
    title('YZ')
end
% -------------------------------------------------------------------------
%% save combined images?
if saveFlag
    imwrite(im_xy_final, fullfile(savePath, 'combined_im_XY.png'))
    imwrite(im_xz_final, fullfile(savePath, 'combined_im_XZ.png'))
    imwrite(im_yz_final, fullfile(savePath, 'combined_im_YZ.png'))
end