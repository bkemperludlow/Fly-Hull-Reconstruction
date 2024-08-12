function qs_force_struct =quasiSteadyForceAll(wingAngleMat, t, wingSide, ...
    params, bodyCM, bodyYPR, plotFlag, smoothFlag, addedMassFlag)
% -------------------------------------------------------------------------
% function to calculate translational, rotational, and added-mass
% components of aerodynamic forces on flapping wing. Uses the quasi-steady
% model described in:
%   (Dickinson, Lehmann, & Sane, 1999 ; Sane & Dickinson, 2002 ;
%     Whitney & Wood, 2010)
%
% ***NOTE: PREVIOUS VERSIONS THAT I (SCW) HAD WRITTEN USED A STUPID FUCKING
% CONVENTION. HERE WE'RE SHIFTING BACK TO A GOOD ONE (SO **X** DIRECTION IS
% ALONG LONG BODY AXIS, RATHER THAN Y)
%
% INPUTS:
%   -wingAngMat: Nx3 matrix containing wing Euler angles in the body frame,
%   in the order [stroke, deviation, (wing) pitch]. *SHOULD BE IN RADIANS
%
%   -t: Nx1 vector with time of each frame (in seconds)
%
%   -wingSide: char indicating whether we're doing the right or left wing
%   (either 'R' or 'L')
%
%   -params: stucture containing wing morphological parameters. This
%   structure can be set to default values via XX.m, but can be adjusted as
%   necessary to accomodate different wing types
%
%   -bodyCM: Nx3 matrix contaning the body position in lab frame. Should be
%   in meters / second
%
%   -bodyYPR: Nx3 matrix containing body Euler angles in the order [yaw,
%   pitch, roll].

% OUTPUT:
%   -qs_force_struct: data structure containing force data + other useful
%   info
%
% -------------------------------------------------------------------------
%% inputs and params
if ~exist('params','var') || isempty(params)
    params = defineQuasiSteadyParams() ;
end
if ~exist('bodyCM','var') || isempty(bodyCM)
    bodyCM = [] ;
end
if ~exist('bodyYPR','var') || isempty(bodyYPR)
    bodyYPR = [] ;
end
if ~exist('plotFlag','var') || isempty(plotFlag)
    plotFlag = false ;
end
if ~exist('smoothFlag','var') || isempty(smoothFlag)
    smoothFlag = false ;
end
if ~exist('addedMassFlag','var') || isempty(addedMassFlag)
    addedMassFlag = false ;
end

% which method to use for calculating angle of attack (alpha)
alphaMode = 'bodyFrame' ; % 'bodyFrame' | 'wingFrame'

% fit type for spline (to be used to create fit objects for d/dt)
fitType = params.fitType ; 

% number of points
N_frames = size(wingAngleMat,1) ;

% make sure t is a column vector
if size(t,1) < size(t,2)
    t = t' ;
end

% also make sure wing angle matrix is Nx3
if size(wingAngleMat,1) < size(wingAngleMat,2)
    wingAngleMat = wingAngleMat' ; 
end

% -----------------------------------------------------------------
%% read in params from struct
body_mass = params.body_mass ;
g = params.g ;
DEG2RAD = (pi/180)  ; %degrees to radians conversion

% ---------------------------------------------------------------------
% if we're calculating for the right wing, need to swap the sign on some
% stuff
if contains(wingSide,'R','IgnoreCase',true)
    if (nanmean(sign(wingAngleMat(:,1))) > 0)
        wingSign = -1 ;
    else
        wingSign = 1 ;
    end
elseif contains(wingSide,'L','IgnoreCase',true)
    wingSign = +1 ;
else
    fprintf('Invalid wingSide entry: %s \n', wingSide)
    keyboard
end

% flip sign of stroke angle (phi) if need be
wingAngleMat(:,1) = wingSign*wingAngleMat(:,1) ;

% make sure wing angles are in radians
strokeMax = nanmedian(abs(wingAngleMat(:,1))) ; 
if (strokeMax > 2*pi)
   wingAngleMat = DEG2RAD.*wingAngleMat ; 
end

% ---------------------------------------------------
%% get span and chord hat vectors (body coordinates)
% if left wing, need to correct sign of rotation angle
wingAngleMatMod = wingAngleMat ; 
if contains(wingSide,'L','IgnoreCase',true)
    wingAngleMatMod(:,3) = pi -  wingAngleMatMod(:,3) ; 
end
[spanHat, chordHat] = getWingVecsFromAng(wingAngleMatMod) ;

% check apan and chord?
if (0)
     h_vec1 = plotWingVectors(1:5:50, spanHat, chordHat, params) ; 
     h_vec2 = plotWingVectors(51:5:100, spanHat, chordHat, params) ; 
end

% ---------------------------------------------------------------
%% calculate wing tip velocity 
% (in lab frame, but expressed in body frame coordinates)
[U_t, omegaWingBody, ~, ~, omegaWingWing] = ...
    calcWingTipVelocity(wingAngleMat, t, params, wingSide, bodyCM, bodyYPR);

% check wing frame angular velocity?
if (0)
    figure
    subplot(2,1,1)
   plot(t, omegaWingWing) ;  
   xlabel('t (s)')
   ylabel('\omega_{ww}')
   subplot(2,1,2)
   plot(t, myNorm(omegaWingWing),'k-') ;  
   xlabel('t (s)')
   ylabel('norm \omega_{ww}')
end

% get unit vector for wing tip velocity
Usq = sum(U_t.*U_t,2) ; % squared magnitude
Uhat = U_t ./ sqrt(Usq) ;

% -----------------------------------------------------------
%% get angle of attack
switch alphaMode
    case 'wingFrame'
        % ----------------------------------------------------------------
        % going to calculate it in wing frame, where it's  atan2(vz, vy)
        
        % pull out wing Euler angles
        phi = wingAngleMat(:,1) ;
        theta = wingAngleMat(:,2) ;
        psi = wingAngleMat(:,3) ;
        
        % write out sines and cosines to make things more compact AND 
        % switch sign on theta (to match Goldstein)
        cph=cos(phi)    ; sph=sin(phi)   ;
        cth=cos(-theta) ; sth=sin(-theta);
        cps=cos(psi)    ; sps=sin(psi)   ;
        
        % get body to wing transformations
        body2wing_all = zeros(3,3,N_frames) ;
        U_t_wing = zeros(size(U_t)) ;
%         testSpan = spanHat ; 
%         testChord = chordHat ; 
        
        % loop over frames
        for i = 1:N_frames
            % generate coordinate transform matrix
            body2wing = [cth(i)*cph(i) ,     cth(i)*sph(i) ,  (-sth(i)) ; ...
                (sps(i)*sth(i)*cph(i)-cps(i)*sph(i)), ...
                (sps(i)*sth(i)*sph(i)+cps(i)*cph(i)),  cth(i)*sps(i) ; ...
                (cps(i)*sth(i)*cph(i)+sps(i)*sph(i)), ...
                (cps(i)*sth(i)*sph(i)-sps(i)*cph(i)), cth(i)*cps(i) ] ;
            

            % write tip velocity in wing frame coordinates
            U_t_wing(i,:) = body2wing*U_t(i,:)' ;
            
%             % test body to wing transformation with span/chord hat
%             testSpan(i,:) = body2wing*testSpan(i,:)' ; 
%             testChord(i,:) = body2wing*testChord(i,:)' ; 
            % store rotation matrix in case we need it later
            body2wing_all(:,:,i) = body2wing ;
        end
        
        % again, applying rotation matrices introduces noise--smooth
        if smoothFlag
            for dim = 1:3
                U_t_wing(:,dim) = smooth(U_t_wing(:,dim),3) ;
            end
        end
        %alpha = atan2(U_t_wing(:,3), U_t_wing(:,2)) ;
        %alpha = atan2(-1.*U_t_wing(:,2), U_t_wing(:,3)) ;
        alpha = atan2(sign(U_t_wing(:,3)).*U_t_wing(:,3), U_t_wing(:,2) ) ;
       
        
    case 'bodyFrame'
        % calculate alpha in body frame by taking angle between velocity
        % vector (relative wind) and chord
        alpha = acos(dot(Uhat, chordHat,2)) ;
        
    otherwise
        fprintf('Invalid mode for angle of attack: %s \n',alphaMode)
        keyboard
end

if smoothFlag
    alpha = smooth(alpha,3) ;
end

% -------------------------------------------------------------------
%% calculate components then sum to get total force
[F_tot, F_t, F_rot, F_a, Lhat, Dhat, wingNormalHat] = ...
    calcQuasiSteadyTerms(U_t, spanHat, chordHat, omegaWingWing, ...
      omegaWingBody, alpha, params, addedMassFlag, smoothFlag) ;

% -------------------------------------------------------------------
%% store data in structure
qs_force_struct = struct() ;

% forces
qs_force_struct.F_tot = F_tot ;                 % total force
qs_force_struct.F_trans = F_t ;                 % translational force component
qs_force_struct.F_rot = F_rot ;                 % rotational force component
qs_force_struct.F_a = F_a ;                     % added mass component

% vectors
qs_force_struct.Lhat = Lhat ;                   % direction of trans. lift force
qs_force_struct.Dhat = Dhat ;                   % direction of trans. drag force
qs_force_struct.Uhat = Uhat ;                   % direction of wingtip velocity
qs_force_struct.wingNormalHat = wingNormalHat ; % wing normal vector
qs_force_struct.spanHat = spanHat ;             % wing span vector
qs_force_struct.chordHat = chordHat ;           % wing chord vector

% angular velocities
qs_force_struct.omegaWingBody = omegaWingBody ; % angular velocity of wing in body frame
qs_force_struct.omegaWingWing = omegaWingWing ; % angular velocity of wing in wing frame
% scalar variables
qs_force_struct.alpha = alpha ;                 % angle of attack
qs_force_struct.t = t ;                         % time (seconds)

% body coordinates
qs_force_struct.bodyCM = bodyCM ;               % body center of mass position
qs_force_struct.bodyYPR = bodyYPR ;             % body euler angle coordinates

% params
qs_force_struct.params = params ;               % general parameters
qs_force_struct.wingSide = wingSide ;           % L or R wing
qs_force_struct.alphaMode = alphaMode ;         % method used to find angle of attack
qs_force_struct.smoothFlag = smoothFlag ;       % whether or not we smoothed angles
qs_force_struct.fitType = fitType  ;            % fit type used for calculating derivatives

% -------------------------------------------------------------------
%% plot results?
if plotFlag
    figure ;
    hold on
    
    % total force
    h_tot = plot(t, myNorm(F_tot)./(body_mass*g), 'k-','LineWidth',1.5) ;
    
    % translational force
    h_trans = plot(t, myNorm(F_t)./(body_mass*g),'r--') ;
    
    % rotational force
    h_rot = plot(t, myNorm(F_rot)./(body_mass*g),'b--') ;
    
    % added mass force
    h_a = plot(t, myNorm(F_a)./(body_mass*g), 'm--') ;

    ylabel('F/mg')
    xlabel('Time (s)')
    
    legend([h_tot, h_trans, h_rot, h_a], ...
        {'Total', 'Trans.', 'Rot.', '+mass'})
    % keyboard
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HELPER FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
%% function to plot wing vectors in 3D as sanity check
function h_main = plotWingVectors(ind, spanHat, chordHat, params ,...
    U_t, LiftVec, DragVec)
% ----------------------------
% inputs
if ~exist('U_t','var') || isempty(U_t)
    U_t = nan(size(spanHat)) ;
end
if ~exist('LiftVec','var') || isempty(LiftVec)
    LiftVec = nan(size(spanHat)) ;
end
if ~exist('DragVec','var') || isempty(DragVec)
    DragVec = nan(size(spanHat)) ;
end


span = params.span ;
chord = params.chord ;
scale = 3e-4 ;
scale2 = 1e-3 ;

N_pts = length(ind) ;

h_main = figure ;
hold on

grays_vec = linspace(0, 0.8, N_pts) ;
grays_vec = fliplr(grays_vec) ;
grays_mat = repmat(grays_vec',1,3) ;

for k = 1:N_pts
    ind_curr = ind(k) ;
    
    % wing span
    h_s = plot3([0.5*span*spanHat(ind_curr,1), span*spanHat(ind_curr,1)],...
        [0.5*span*spanHat(ind_curr,2), span*spanHat(ind_curr,2)],...
        [0.5*span*spanHat(ind_curr,3), span*spanHat(ind_curr,3)],...
        '-','LineWidth',2,'Color',grays_mat(k,:)) ;
    
    % wing chord
    h_c = plot3(0.75*span*spanHat(ind_curr,1) + [0, chord*chordHat(ind_curr,1)],...
        0.75*span*spanHat(ind_curr,2) + [0, chord*chordHat(ind_curr,2)],...
        0.75*span*spanHat(ind_curr,3) + [0, chord*chordHat(ind_curr,3)],...
        'r-','LineWidth',2) ;
    
    % wing velocity
    h_v = plot3(span*spanHat(ind_curr,1) + [0, scale*U_t(ind_curr,1)],...
        span*spanHat(ind_curr,2) + [0, scale*U_t(ind_curr,2)],...
        span*spanHat(ind_curr,3) + [0, scale*U_t(ind_curr,3)],...
        'b-','LineWidth',2) ;
    % wing drag
    h_d = plot3(0.75*span*spanHat(ind_curr,1) + [0, scale2*DragVec(ind_curr,1)],...
        0.75*span*spanHat(ind_curr,2) + [0, scale2*DragVec(ind_curr,2)],...
        0.75*span*spanHat(ind_curr,3) + [0, scale2*DragVec(ind_curr,3)],...
        'c-','LineWidth',2) ;
    
    
    % wing lift
    h_l = plot3(0.75*span*spanHat(ind_curr,1) + [0, scale2*LiftVec(ind_curr,1)],...
        0.75*span*spanHat(ind_curr,2) + [0, scale2*LiftVec(ind_curr,2)],...
        0.75*span*spanHat(ind_curr,3) + [0, scale2*LiftVec(ind_curr,3)],...
        'm-','LineWidth',2) ;
    
end

axis equal
grid on
box on

xlabel('X')
ylabel('Y')
zlabel('Z')

legend([h_s, h_c, h_v, h_d, h_l], {'Span', 'Chord', 'Tip Velocity',...
    'Drag','Lift'}) ;

end


% % -------------------------------------------------------------------------
% %% get angle of attack (from Tsevi's code)
% function [alpha, U] = getAngleOfAttack(phi, theta, psi, phidot, thetadot, ...
%     psidot, vBody, vCent, params)
%
% % --------------------------------------
% % get wing morphological parameters
% span = params.span ;
% chord = params.chord ;
%
% % -------------------------------------------
% % wing angle sines and cosines
% cph=cos(phi)    ; sph=sin(phi)   ;
% cth=cos(-theta) ; sth=sin(-theta);
% cps=cos(psi)    ; sps=sin(psi)   ;
%
% % transformation between body and wing frames
% body2wing = [cth*cph            cth*sph           (-sth) ; ...
%     (sps*sth*cph-cps*sph) (sps*sth*sph+cps*cph) cth*sps ; ...
%     (cps*sth*cph+sps*sph) (cps*sth*sph-sps*cph) cth*cps ] ;
%
% wing2body = body2wing' ;
%
% % Goldstein appendix = attila's thesis with -1*theta
% %  should have used the name omegaWing rather than omegaBody
% modThetadot=-thetadot ;
% omegaWing = [ psidot - phidot*sth   ; ... % note minus sign is definition of sth above, so minus here is ok
%     modThetadot*cps  + phidot*cth*sps  ; ...
%     -modThetadot*sps + phidot*cth*cps ] ;
%
% % 6-vector velocity in the body frame of reference (bs=body system)
%
% omegaLab = [psidot*cth*cph - modThetadot*sph ; ...
%     psidot*cth*sph  + modThetadot*cph ; ...
%     phidot - psidot*sth] ; % XXX
%
%
% % -------------------------------------
% % body velocity
% if (isempty(vBody))
%     vBody = 0 ;
% end
%
% % --------------------------------------------------
% % translational velocity of wing CM
% if (isempty(vCent))
%     % http://planning.cs.uiuc.edu/node102.html
%     % rotation matrices
%     %{
%     R1 = [ cph -sph 0 ; ...
%         sph cph 0 ; ...
%         0    0  1 ] ;
%
%     R2 = [cth 0 sth ; ...
%         0     1  0 ; ...
%         -sth  0  cth] ;    %
%     %}
%     rcent = (wing2body)*[span/2 ;-chord/2 ;0] ;     % wing center position
%     vlab = cross(omegaLab, rcent) + vBody; % velocity of wing center in the lab frame
% else
%     vlab = vCent ;
% end
%
% vWingFrame = body2wing*vlab ; % wing velocity in the wing frame (vBody is included here)
%
% vbs = [vWingFrame ; omegaWing ];
%
% % --------------------------------------------------------------
% % use wing translational and angular velocity to get alpha
% rv    = zeros(4, 1) ;
% rv(1) = span / 2; % previously was span/2.
% rv(4) = 1 ;
%
% v = wedge(vbs) * rv ;  % v = cross(omega, span) + translational velocity
%
% vyp = v(2) ; % velocity component perpendicular to wing surface
% vzp = v(3) ; % velocity component parallel to chord vector
%
% U     = [ 0, vyp, vzp] ;
%
% alpha = atan2(vzp, vyp) ;
%
%
% end
%

%{
% visually check direction of normal vector
if (0)
   F_rot_hat = F_rot./myNorm(F_rot) ; 
   figure ; 
    subplot(2,1,1)
    hold on
    plot3(spanHat(1:50,1), spanHat(1:50,2), spanHat(1:50,3),'k.-') ; 
    plot3(spanHat(1,1), spanHat(1,2), spanHat(1,3),'ko','markerfacecolor',0.7*[0,1,0])
    for ii = 1:2:50 
        plot3(spanHat(ii,1) + 0.25*[0, chordHat(ii,1)],...
            spanHat(ii,2) + 0.25*[0, chordHat(ii,2)], ...
            spanHat(ii,3) + 0.25*[0, chordHat(ii,3)],'b.-') ; 
    end
    for ii = 1:2:50  
        plot3(spanHat(ii,1) + 0.25*[0, F_rot_hat(ii,1)], ...
            spanHat(ii,2) + 0.25*[0, F_rot_hat(ii,2)], ...
            spanHat(ii,3) + 0.25*[0, F_rot_hat(ii,3)],'ro-') ; 
    end
    axis equal ; grid on
    title('Fwd')
    
    subplot(2,1,2)
    hold on
    plot3(spanHat(51:end,1), spanHat(51:end,2), spanHat(51:end,3),'k.-') ; 
    plot3(spanHat(51,1), spanHat(51,2), spanHat(51,3),'ko','markerfacecolor',0.7*[0,1,0])
    for ii = 51:2:N_frames
        plot3(spanHat(ii,1) + 0.25*[0, chordHat(ii,1)],...
            spanHat(ii,2) + 0.25*[0, chordHat(ii,2)], ...
            spanHat(ii,3) + 0.25*[0, chordHat(ii,3)],'b.-') ; 
    end
    for ii = 51:2:N_frames
        plot3(spanHat(ii,1) + 0.25*[0, F_rot_hat(ii,1)], ...
            spanHat(ii,2) + 0.25*[0, F_rot_hat(ii,2)], ...
            spanHat(ii,3) + 0.25*[0, F_rot_hat(ii,3)],'ro-') ; 
    end
    axis equal ; grid on
    title('Back')
end
%}