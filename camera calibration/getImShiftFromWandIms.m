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
function [centers_3D, im_shifts, centers_2D] = getImShiftFromWandIms(wandImageStruct,...
    dlt, easyWandData, voxelSize, volLength, order, debugFlag)
%-----------------------
%% params and inputs
if ~exist('voxelSize','var') || isempty(voxelSize)
    voxelSize = 5.0e-5 ;
end
if ~exist('volLength','var') || isempty(volLength)
    volLength = 0.008 ;
end
if ~exist('order','var') || isempty(order)
    order = [2, 1, 3] ;
end
if ~exist('debugFlag','var') || isempty(debugFlag)
    debugFlag = false ;
end

% which method to use to solve for translations. having trouble with non
% smooth problems :(
solveType = 'fmincon' ; % 'brute_force', 'fmincon', 'fminunc','fminsearch', 'particle_swarm'

YZ = 1 ;
XZ = 2 ;
XY = 3 ;

sphere_types = {'big', 'small'} ;
cam_names = {'yz', 'xz', 'xy'} ;

imSize =  512 ;
N_images = length(wandImageStruct) ;
N_cams = length(cam_names) ;
N_spheres = length(sphere_types) ;
% --------------------------------
%% initialize data structure
centers_3D = nan(N_images, 3, N_spheres) ;
im_shifts = nan(N_images, 2, N_spheres, N_cams-1) ;
%centers_2D = nan(N_images, N_spheres, N_cams, 2) ; 
% ------------------------------------------
%% loop through images/sphere types/cameras
% ------------------------------------------
for i = 1:N_images
    for j = 1:N_spheres
        sphere_type = sphere_types{j} ;
        
        % -------------------------------------------
        %% get center of mass and radius for sphere in all views
        centers_curr = nan(N_cams, 2) ;
        radii_curr = nan(N_cams, 1) ; 
        for k = 1:N_cams
            % center
            center_temp = ...
                wandImageStruct(i).([cam_names{k} '_center_' sphere_type]) ;
            centers_curr(k,:) = [center_temp(1), imSize - center_temp(2)] ;
            
            % radius
            radius_temp = ...
                wandImageStruct(i).([cam_names{k} '_radius_' sphere_type]) ;
            radii_curr(k) = radius_temp ; 
        end
        %centers_2D(i,j,:,:) = centers_curr ; 
        
        % get 3D center point
        Ctr = getFlyCenterSeed(centers_curr(1,:), centers_curr(2,:), ...
            centers_curr(3,:), easyWandData, dlt, order, voxelSize) ;
        
        centers_3D(i,:,j) = Ctr ; 
        % -------------------------------------
        %% generate voxel space for this frame
        % find the [X0, Y0, Z0] which is the real-space coordinates of the corner
        %  of the current sampling sub-volume.
        X0 = Ctr(1) - volLength/2 ; % assuming that volLength/2 is an integer number of voxelSize
        Y0 = Ctr(2) - volLength/2 ;
        Z0 = Ctr(3) - volLength/2 ;
        
        % generate the real-space voxel grid vectors X, Y, Z
        % (previous code started from the middle but the X,Y,Z were the same).
        X = X0 : voxelSize : X0+volLength ;
        Y = Y0 : voxelSize : Y0+volLength ;
        Z = Z0 : voxelSize : Z0+volLength ;
        [Xgrid, Ygrid, Zgrid] = meshgrid(X, Y, Z) ;
        points = [Xgrid(:) , Ygrid(:), Zgrid(:)] ;
        clear Xgrid Ygrid Zgrid
        % find the offset of the real space volume from the real space origin
        % expressed in no. of voxels
        offset = int16([X0, Y0, Z0]/voxelSize); % +++
        
        if ( double(offset(1))*voxelSize - X0 > eps || ...
                double(offset(2))*voxelSize - Y0 > eps || ...
                double(offset(3))*voxelSize - Z0 > eps)
            error('something is wrong with voxel or offset calcuation') ;
        end
        
        % size of voxel space
        Nx = length(X) ; Ny = length(Y) ; Nz = length(Z) ;
        voxSizeVec = [Nx, Ny, Nz] ;
        
        % --------------------------------------
        %% fit for best image translation
        % first get xy projection
        im_xy = drawBWcircle(centers_curr(XY,:), radii_curr(XY)) ; 
        isXY = checkVoxelCamSphere(points, dlt(:,order(XY)), im_xy, imSize) ;
        % also get side view images and dlt coeffs
        im_yz = drawBWcircle(centers_curr(YZ,:), radii_curr(YZ)) ; 
        im_xz = drawBWcircle(centers_curr(XZ,:), radii_curr(XZ)) ; 
        dlt_yz = dlt(:, order(YZ)) ;
        dlt_xz = dlt(:, order(XZ)) ;
        
        % set initial guesses for solver. x is a vector of image
        % translations with:
        % x = [shift_yz_horz, shift_yz_vert, shift_xz_horz, shift_xz_vert]
        x0 = [0; 0; 0; 0] ;
        
        % --------------------------
        % run optimization routine
        switch solveType
            case 'brute_force'
                % nb: for now, only looking at horizontal translations to
                % save (some) time
                x1 = (x0(1)-5) : (x0(1)+5) ;
                x2 = 0 ;
                x3 = (x0(3)-5) : (x0(3)+5) ;
                x4 = 0 ;
                fval_mat = nan(length(x1), length(x3), length(x4)) ;
                cc = 1 ;
                for ii = 1:length(x1)
                    for jj = 1:length(x3)
                        for kk = 1:length(x4)
                            x = [x1(ii), x2, x3(jj), x4(kk)] ;
                            fval_mat(ii, jj, kk) = wandShiftObjective(x,...
                                im_yz, im_xz, dlt_yz, dlt_xz, points, ...
                                isXY, imSize) ;
                            disp(cc)
                            cc = cc + 1 ;
                        end
                    end
                end
                
                [~, min_idx] = min(fval_mat(:)) ;
                [a, b, c] = ind2sub(size(fval_mat), min_idx) ;
                x = [x1(a), x2, x3(b), x4(c)] ;
            case 'fmincon'
                A = eye(length(x0)) ;
                b = 10*ones(1,length(x0)) ;
                opt = optimoptions(@fmincon,'Display','iter',...
                    'OptimalityTolerance',1e-10,...
                    'FiniteDifferenceType','central') ;
                [x, fval] = fmincon(@(x)wandShiftObjective(x,im_yz, im_xz,...
                    dlt_yz,dlt_xz, points, isXY, imSize), x0, A, b, [],[],...
                    [],[],[],opt) ;
            case 'fminunc'
                opt = optimoptions(@fminunc,'Display','iter',...
                    'OptimalityTolerance',1e-10) ; 
                [x, fval] = fminunc(@(x)wandShiftObjective(x,im_yz, ...
                    im_xz, dlt_yz, dlt_xz, points, isXY, imSize), x0, opt) ;
            case 'fminsearch'
                opt = optimoptions(@fminsearch,'Display','iter') ; 
                [x, fval] = fminsearch(@(x)wandShiftObjective(x,im_yz, ...
                    im_xz, dlt_yz, dlt_xz, points, isXY, imSize), x0, opt) ;
            case 'particle_swarm'
                x = particleswarm(@(x)wandShiftObjective(x,im_yz, im_xz, ...
                    dlt_yz, dlt_xz, points, isXY, imSize),4) ; 
            otherwise
                disp('invalid solver selection')
                keyboard 
        end

        % convert minimizer output into image translations
        im_shift_yz = [x(1), x(2)] ;
        im_shift_xz = [x(3), x(4)] ;
        
        % --------------------------
        %% check results?
        if debugFlag
            % voxel plot
            h_vox = figure ;
            
            % find shifted and old hulls
            im_yz_shifted = imtranslate(im_yz, im_shift_yz) ;
            im_xz_shifted = imtranslate(im_xz, im_shift_xz) ;
            
            isYZ_shift = checkVoxelCamSphere(points, dlt_yz, ...
                im_yz_shifted, imSize) ;
            isXZ_shift = checkVoxelCamSphere(points, dlt_xz, ...
                im_xz_shifted, imSize) ;
            
            isYZ_old = checkVoxelCamSphere(points, dlt_yz, ...
                im_yz, imSize) ;
            isXZ_old = checkVoxelCamSphere(points, dlt_xz, ...
                im_xz, imSize) ;
            
            % find voxels seen by all 3 views
            hull_idx_shift = all([isXY, isYZ_shift, isXZ_shift],2) ;
            hull_idx = all([isXY, isYZ_old, isXZ_old],2) ;
            
            % --------------
            subplot(1,2,1)
            plot3(points(hull_idx,1), points(hull_idx,2), ...
                points(hull_idx,3), 'rx')
            axis equal
            box on
            grid on
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            
            % --------------
            subplot(1,2,2)
            plot3(points(hull_idx_shift,1), points(hull_idx_shift,2), ...
                points(hull_idx_shift,3), 'b.')
            
            axis equal
            box on
            grid on
            xlabel('X')
            ylabel('Y')
            zlabel('Z')
            % ---------------------
            % image plot
            h_im = figure ;
            subplot(1,2,1)
            imshowpair(im_yz, im_yz_shifted)
            title('YZ')
            
            subplot(1,2,2)
            imshowpair(im_xz, im_xz_shifted)
            title('XZ')
        end
        
        if any(abs(x) > eps )
            keyboard
        end
        %% store results
        centers_3D(i, :, j) = Ctr ;
        im_shifts(i, :, j, 1) = im_shift_yz;
        im_shifts(i, :, j, 2) = im_shift_xz;
        
        
    end
    fprintf('Completed frame %d /%d \n',i, N_images)
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
function Ctr = getFlyCenterSeed(center_mass_yz, center_mass_xz, ...
    center_mass_xy, easyWandData, dlt, order, voxelSize)
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

%R=change_ref_frame(Rxz,Rxy);
%Ctr=R*Ctr;
% ---------------------------------------------
% 3D HULL RECONSTRUCTION - FIND THE BODY VOXELS
% ---------------------------------------------


% shift Ctr to be the middle of some voxel.
% given that volLength/2 / params.voxelSize is integer
% NOTE: the points (0,0,0) in real space is the center of the camera
% system, i.e. the center-of-mass of all the points used for
% calibration.

Ctr(1) = round(Ctr(1)/voxelSize) * voxelSize ;
Ctr(2) = round(Ctr(2)/voxelSize) * voxelSize ;
Ctr(3) = round(Ctr(3)/voxelSize) * voxelSize ;
end
% -------------------------------------------------------------------------
%% check voxel against image
function isVoxOn = checkVoxelCamSphere(r, dlt, imSphere, detectorLengthPix)
% project all pixels using dlt, then round
pixel   = dlt_inverse(dlt, r );
pixel(:,1) = round(pixel(:,1)) ;
pixel(:,2) = round(pixel(:,2));

% take only pixels that are within the range of the image (i.e. 0 < x <
% imageWidth)
good_idx = (pixel(:,1) > 0) & (pixel(:,1) < detectorLengthPix) & ...
    (pixel(:,2) > 0) & (pixel(:,2) < detectorLengthPix) ;
% convert these to valid entries for now, but we'll switch them later
pixel(~good_idx,:) = 1 ;

% get linear index for pixels to test
pixel_idx = sub2ind(size(imSphere), pixel(:,2), pixel(:,1)) ;

% test pixels against images
isVoxOn = imSphere(pixel_idx) ;
isVoxOn(~good_idx) = false ;
end
% -------------------------------------------------------------------------
%% wand shift objective function
function [f,g,H] = wandShiftObjective(xin, im_yz, im_xz, dlt_yz, dlt_xz,...
    points, isXY, imSize)

% read input vector into image translations
im_shift_yz = [xin(1), xin(2)] ;
im_shift_xz = [xin(3), xin(4)] ;

% translate image
im_yz_shifted = imtranslate(im_yz, im_shift_yz) ;
im_xz_shifted = imtranslate(im_xz, im_shift_xz) ;

% find the indices of points projecting onto yz and xz
isYZ = checkVoxelCamSphere(points, dlt_yz, im_yz_shifted, imSize) ;
isXZ = checkVoxelCamSphere(points, dlt_xz, im_xz_shifted, imSize) ;

% find voxels seen by all 3 views, get size
hull_idx = all([isXY, isYZ, isXZ],2) ;
hull_size = sum(hull_idx) ;

% calculate output of objective function (need to make it a minimization
% thing)
f = 1/(1+hull_size) ;
%f = -1*hull_size ;

% jacobian and hessian; probably won't ever do
if nargout > 1
    g = zeros(size(x)) ;
    if nargout > 2
        H = zeros(length(x),length(x)) ;
    end
end
end
