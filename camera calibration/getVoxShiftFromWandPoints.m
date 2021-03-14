%--------------------------------------------------------------------------
% attempt to use the wand images to find image transformations that will
% refine reconstruction.
%
% the basic idea is that, for each frame where the sphere is visible in all
% views, we should be able to find the image transformations for each of
% the camera views that maximizes the hull size of the reconstructed
% sphere. this will give us a sampling of 3D space that we can interpolate
% through and calculate shifts for any positon that the fly is in
%--------------------------------------------------------------------------
function [centers_3D, vox_shifts, centers_2D] = ...
    getVoxShiftFromWandPoints(calibrationPath, voxelSize,...
        order, debugFlag, saveFlag)
%-----------------------
%% params and inputs
if ~exist('voxelSize','var') || isempty(voxelSize)
    voxelSize = 5.0e-5 ;
end
if ~exist('order','var') || isempty(order)
    order = [2, 1, 3] ;
end
if ~exist('debugFlag','var') || isempty(debugFlag)
    debugFlag = true ;
end
if ~exist('saveFlag','var') || isempty(saveFlag)
    saveFlag = true ;
end

% load in easyWandData and calibration coefficients
DLT_filename = fullfile(calibrationPath, 'calibration_dltCoefs.csv') ;
easyWand_filename = fullfile(calibrationPath, 'calibration_easyWandData.mat') ;
dlt = load(DLT_filename) ; % CSV
load(easyWand_filename); % contains easyWandData
set(0, 'ShowHiddenHandles', 'on')
close 'easyWand 5'
wandPoints = easyWandData.wandPts ; 
% camera indices
YZ = 1 ;
XZ = 2 ;
XY = 3 ;

sphere_types = {'big', 'small'} ;
cam_names = {'yz', 'xz', 'xy'} ;

imSize = 512 ;
N_images = size(wandPoints,1) ;
N_cams = length(cam_names) ;
N_spheres = length(sphere_types) ;
% --------------------------------
%% initialize data structure
centers_3D = nan(N_images, 3, N_spheres) ;
vox_shifts = nan(N_images, 3, N_spheres, N_cams-1) ;
centers_2D = nan(N_images, N_spheres, N_cams, 2) ; 

% ------------------------------------------
%% loop through images/sphere types/cameras
% ------------------------------------------
for i = 1:N_images
    for j = 1:N_spheres
        % -------------------------------------------
        %% get center of mass for sphere in all views
        centers_curr = nan(N_cams, 2) ;
        %radii_curr = nan(N_cams, 1) ; 
        for k = 1:N_cams
            % center
            idx1 = (2*order(k) - 1) + 6*mod(j,2) ; 
            idx2 = 2*order(k) + 6*mod(j,2) ; 
            center_temp = wandPoints(i, idx1:idx2) ; 
            centers_curr(k,:) = [center_temp(1), center_temp(2)] ;
        end
        centers_2D(i,j,:,:) = centers_curr ; 
        
        % --------------------------------------
        %% get 3D center point and voxel shifts
        [Ctr, YZ_shift, XZ_shift] = getCenterAndShift(centers_curr(1,:),...
            centers_curr(2,:), centers_curr(3,:), easyWandData, dlt,...
            order, voxelSize) ;
        
        centers_3D(i,:,j) = Ctr ; 
        
        vox_shifts(i, :, j, YZ) = YZ_shift ;
        vox_shifts(i, :, j, XZ) = XZ_shift ;
        
        
    end
    %fprintf('Completed frame %d /%d \n',i, N_images)
end

% -------------------------------------------------------------------------
%% just concatenate results--not sure why i separated in the first place
centers_3D = vertcat(squeeze(centers_3D(:,:,1)), squeeze(centers_3D(:,:,2))) ;
vox_shifts = vertcat(squeeze(vox_shifts(:,:,1,:)), squeeze(vox_shifts(:,:,2,:))) ;
% -------------------------------------------------------------------------
%% make plot to show results
if debugFlag 
    h_shifts = figure('PaperPositionMode','auto') ; 
    for i = 1:2
       
       vox_shifts_curr = squeeze(vox_shifts(:,:,i)) ;
       shifts_norm = myNorm(vox_shifts_curr) ; 
       shifts_hat = vox_shifts./repmat(shifts_norm,1,3) ; 
       shifts_ang_el = asin(shifts_hat(:,3)) ; 
       shifts_ang_az = atan2(shifts_hat(:,2), shifts_hat(:,1)) ; 
       
       scatter_size = 1+ 12*(shifts_norm - min(shifts_norm)) ./ ...
           (max(shifts_norm) - min(shifts_norm)) ;
       % -----------
       % az angle
       subplot(2, 2, 1 + (i-1)*2)
       scatter3(centers_3D(:,1), centers_3D(:,2), centers_3D(:,3),...
           scatter_size, shifts_ang_az, 'filled')
       title([ cam_names{i} ' az angle'])
       axis tight
       %axis equal 
       box on 
       grid on
       
       % -----------
       % el angle
       subplot(2, 2, 2 + (i-1)*2)
       scatter3(centers_3D(:,1), centers_3D(:,2), centers_3D(:,3),...
           scatter_size, shifts_ang_el, 'filled')
       title([ cam_names{i} ' el angle'])
       axis tight
       %axis equal 
       box on 
       grid on
    end
end

% -----------------
%% save results?
if saveFlag
   save(fullfile(calibrationPath, 'voxShifts.mat'), 'centers_3D', ...
       'vox_shifts', 'centers_2D')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%% get camera transformation matrices from dlt coefficients
function bw_circle = drawBWcircle(center, radius, imageSizeX, imageSizeY)

if ~exist('imageSizeX','var') || isempty(imageSizeX)
    imageSizeX = 512 ;
end
if ~exist('imageSizeY','var') || isempty(imageSizeY)
    imageSizeY = 512 ;
end

[columnsInImage, rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
bw_circle = (rowsInImage - center(2)).^2 ...
    + (columnsInImage - center(1)).^2 <= radius.^2;
end

%--------------------------------------------------------------------------
%% get camera transformation matrices from dlt coefficients
function [Kxz,Kyz,Kxy,Rxz,Ryz,Rxy] = find_transform_dlt(easyWandData, dlt,order)

%load(calibration_directory)
%dlt=load(coefs_directory);
%order is a 1x3 vector  such that order(1)= index of camera YZ for calibration
%                                 order(2)= index of camera XZ for calibration
%                                 order(3)= index of camera XY for calibration

dltxy=dlt(:,order(3));
dltxz=dlt(:,order(2));
dltyz=dlt(:,order(1));

pp=easyWandData.ppts;
f=easyWandData.focalLengths';

f=f(:,order);
pp=pp(:,[2*order(1)-1,2*order(1), 2*order(2)-1,2*order(2),...
    2*order(3)-1,2*order(3)]);

Kyz=[f(1),0,pp(1);0,f(1),pp(2);0,0,1];
Kxz=[f(2),0,pp(3);0,f(2),pp(4);0,0,1];
Kxy=[f(3),0,pp(5);0,f(3),pp(6);0,0,1];

Axz=[dltxz(1),dltxz(2),dltxz(3),dltxz(4); ...
    dltxz(5),dltxz(6),dltxz(7),dltxz(8); ...
    dltxz(9),dltxz(10),dltxz(11),1];
Ayz=[dltyz(1),dltyz(2),dltyz(3),dltyz(4); ...
    dltyz(5),dltyz(6),dltyz(7),dltyz(8); ...
    dltyz(9),dltyz(10),dltyz(11),1];
Axy=[dltxy(1),dltxy(2),dltxy(3),dltxy(4); ...
    dltxy(5),dltxy(6),dltxy(7),dltxy(8); ...
    dltxy(9),dltxy(10),dltxy(11),1];

Rxz= Kxz \ Axz ; % inv(Kxz)*Axz;
Ryz= Kyz \ Ayz ; % inv(Kyz)*Ayz;
Rxy= Kxy \ Axy ; % inv(Kxy)*Axy;

Rxz=Rxz/max([norm(Rxz(1:3,1)),norm(Rxz(1:3,2)),norm(Rxz(1:3,3))]);
Ryz=Ryz/max([norm(Ryz(1:3,1)),norm(Ryz(1:3,2)),norm(Ryz(1:3,3))]);
Rxy=Rxy/max([norm(Rxy(1:3,1)),norm(Rxy(1:3,2)),norm(Rxy(1:3,3))]);

end
% -------------------------------------------------------------------------
%% get fly center
function [Ctr, YZ_shift, XZ_shift] = getCenterAndShift(center_mass_yz, ...
    center_mass_xz, center_mass_xy, easyWandData, dlt, order, voxelSize)
%---------------------------------------------------------------------
% first get camera transformation matrices
[Kxz,Kyz,Kxy,Rxz,Ryz,Rxy] = find_transform_dlt(easyWandData, dlt,order) ;
%---------------------------------------------------------------------
% For each pair of camera/image find the camera center 3d coordinates
% and a direction vector of the 3d line that connects
% the camera center i with the center of mass from the image of camera i
%-----------------------------------------------------------------------
Uyz =   Ryz(1:3,1:3)' * (Kyz \ ([center_mass_yz,1]')); % Ryz(1:3,1:3)' * inv(Kyz) * [center_mass_yz,1]';
Oyz = - Ryz(1:3,1:3)' * Ryz(1:3,4);

Uxz =   Rxz(1:3,1:3)' * (Kxz \ ([center_mass_xz,1]')); % Rxz(1:3,1:3)' * inv(Kxz) * [center_mass_xz,1]';
Oxz = - Rxz(1:3,1:3)' * Rxz(1:3,4);

Uxy =   Rxy(1:3,1:3)' * (Kxy \ ([center_mass_xy,1]')); % Rxy(1:3,1:3)' * inv(Kxy) * [center_mass_xy,1]';
Oxy = - Rxy(1:3,1:3)' * Rxy(1:3,4);

% Find G1 such that d(G1,line xz) and d(G1,line yz) are minimal

A1=[-dot(Uyz,Uyz),dot(Uyz,Uxz);-dot(Uyz,Uxz),dot(Uxz,Uxz)];
B1=[-dot(Oxz-Oyz,Uyz);-dot(Oxz-Oyz,Uxz)];

lambda1=linsolve(A1,B1);

G1= (Oyz+Oxz+lambda1(1)*Uyz+lambda1(2)*Uxz)/2;

% Find G2 such that d(G1,line xz) and d(G1,line xy) are minimal

A2=[-dot(Uxy,Uxy),dot(Uxy,Uxz);-dot(Uxy,Uxz),dot(Uxz,Uxz)];
B2=[-dot(Oxz-Oxy,Uxy);-dot(Oxz-Oxy,Uxz)];

lambda2=linsolve(A2,B2);

G2= (Oxy+Oxz+lambda2(1)*Uxy+lambda2(2)*Uxz)/2;

% Find G3 such that d(G3,line yz) and d(G3,line xy) are minimal

A3 = [-dot(Uxy,Uxy),dot(Uxy,Uyz);-dot(Uxy,Uyz),dot(Uyz,Uyz)];
B3 = [-dot(Oyz-Oxy,Uxy);-dot(Oyz-Oxy,Uyz)];

lambda3 = linsolve(A3,B3);

G3 = (Oxy+Oyz+lambda3(1)*Uxy+lambda3(2)*Uyz)/2;

% Take as center the mean of the three points.
G = [G1';G2';G3'];

Ctr = mean(G); % estimated c.m. in real space.

% shift Ctr to be the middle of some voxel.
% given that volLength/2 / params.voxelSize is integer
% NOTE: the points (0,0,0) in real space is the center of the camera
% system, i.e. the center-of-mass of all the points used for
% calibration.

Ctr(1) = round(Ctr(1)/voxelSize) * voxelSize ;
Ctr(2) = round(Ctr(2)/voxelSize) * voxelSize ;
Ctr(3) = round(Ctr(3)/voxelSize) * voxelSize ;

% ------------------------------------------------------
% get shifts to make centers align as well as possible
[c1_yz, c2_yz] = shortestConnLine(Uyz, Oyz, Uxy, Oxy) ; 
YZ_shift = c2_yz - c1_yz ; 
[c1_xz, c2_xz] = shortestConnLine(Uxz, Oxz, Uxy, Oxy) ; 
XZ_shift = (c2_xz - c1_xz) + (c2_yz - c2_xz) ; 

if (0)
   figure ; 
   hold on
   plot3([Oyz(1), Uyz(1)], [Oyz(2), Uyz(2)], [Oyz(3), Uyz(3)],'c--','lineWidth',2)
   plot3([Oxz(1), Uxz(1)], [Oxz(2), Uxz(2)], [Oxz(3), Uxz(3)],'m--','lineWidth',2)
   plot3([Oxy(1), Uxy(1)], [Oxy(2), Uxy(2)], [Oxy(3), Uxy(3)],'b-','lineWidth',2)
   
   plot3([Oyz(1), Uyz(1)] + YZ_shift(1), [Oyz(2), Uyz(2)] + YZ_shift(2), ...
       [Oyz(3), Uyz(3)] + YZ_shift(3),'r-','lineWidth',2)
   plot3([Oxz(1), Uxz(1)] + XZ_shift(1), [Oxz(2), Uxz(2)] + XZ_shift(2), ...
       [Oxz(3), Uxz(3)] + XZ_shift(3),'g-','lineWidth',2)

   box on 
   grid on
   %axis equal
    
end
end
% -------------------------------------------------------------------------
%% shortest line segment joining two skew lines (wikipedia)
% v1, v2 are arbitray points on line v
% u1, u2 are arbitray points on line u
function [c1, c2] = shortestConnLine(v1, v2, u1, u2)

% write lines as vectors
p1 = v1 ; 
d1 = (v2-v1)./norm(v2-v1) ; 

p2 = u1 ; 
d2 = (u2-u1)./norm(u2-u1) ; 

% find new line that is perpendicular to both lines
n = cross(d1, d2) ; 

% plane formed by translations of u along n contains p2 and is perp to n2:
n2 = cross(d2, n) ; 
n1 = cross(d1, cross(d2, d1)) ; 

% find points
c1 = p1 + (dot((p2-p1), n2)/dot(d1,n2)).*d1 ; 
c2 = p2 + (dot((p1-p2), n1)/dot(d2,n1)).*d2 ; 

end
