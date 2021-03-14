function [newChord, newAltChord, newDiag1, newDiag2] = refineChordVector(data, wingLetter)
% find the frames that correspond to forward and back stroke
%
% find the stroke plane
%
% still need to handle cases where the wing points down with respect to the stroke
% plane
%
% also need to handle cases where "diag" difference is very large.
%
SPEED_THRESHOLD = 1; % m/s
DELTA = 17; % take +17 and -17 points from each time point. to find the "current" stroke plane.

plotFlag = true ;

switch upper(wingLetter)
    case 'R'
        % find the wing tip with respect to the body cm
        tip      = (data.rightWingTips - data.bodyCM) * data.params.voxelSize ; % in METERS
        othertip = (data.leftWingTips  - data.bodyCM) * data.params.voxelSize ; % in METERS
        %span = data.rightSpanHats ;
        mainChord = data.rightChordHats ;
        altChord  = data.chord1AltHats ;
        diag1 = data.diag11Right ;
        diag2 = data.diag12Right ;
    case 'L'
        % find the wing tip with respect to the body cm
        tip      = (data.leftWingTips  - data.bodyCM) * data.params.voxelSize ; % in METERS
        othertip = (data.rightWingTips - data.bodyCM) * data.params.voxelSize ; % in METERS
        %span = data.rightSpanHats ;
        mainChord = data.leftChordHats ;
        altChord  = data.chord2AltHats ;
        diag1 = data.diag21Left ;
        diag2 = data.diag22Left ;
    otherwise
        error('wingLetter can be only ''L'' or ''R''.') ;
end

dt      = 1 / data.params.fps ;
N       = size(tip,1) ;
tvec    = (0:N-1)' * dt ;
tvec_ms = tvec * 1000 ;

% find a fit for the stroke plane that goes through "tip" and "othertip"

strokePlaneNormals = zeros(N,3) ;

for it=1:N
    t1 = max([1, it-DELTA]) ;
    if (t1==1)
        t2 = 2*DELTA+1 ;
    else
        t2 = min([N, it+DELTA]);
        if (t2==N)
            t1 = N-2*DELTA ;
        end
    end
    if (t2-t1~=2*DELTA)
        disp('problem') ;
        keyboard ;
    end
    
    [n, V, p] = affine_fit([tip(t1:t2,:) ; othertip(t1:t2,:)]) ;
    strokePlaneNormals(it,:) = n' ;
    
end

% spline smooth the tip position, find tolerance "tol" that gives (+/-) of
% a given ESTERR


ESTERR = data.params.voxelSize * 4; % measurement is (+/-)esterr about the "real" data in meters
tol = N * ESTERR^2 * dt ;

%disp('finding time-points to ignore. this is specific to expr7mov9!!!') ;

[sp_tipx, tipx_smooth] =  spaps(tvec, tip(:,1), tol) ;
[sp_tipy, tipy_smooth] =  spaps(tvec, tip(:,2), tol) ;
[sp_tipz, tipz_smooth] =  spaps(tvec, tip(:,3), tol) ;

sp_tipvx = fnder(sp_tipx, 1) ;
sp_tipvy = fnder(sp_tipy, 1) ;
sp_tipvz = fnder(sp_tipz, 1) ;

tip_vx = fnval(sp_tipvx, tvec) ;
tip_vy = fnval(sp_tipvy, tvec) ;
tip_vz = fnval(sp_tipvz, tvec) ;

tip_v = [tip_vx, tip_vy, tip_vz] ;

% calc wing velocity vector in the current stroke plane
tip_v_plane = tip_v -  strokePlaneNormals .* repmat ( dot(tip_v, strokePlaneNormals, 2) , 1 , 3) ;

speed = (tip_vx.^2 + tip_vy.^2 + tip_vz.^2) .^ 0.5 ;
speed_plane = myNorm(tip_v_plane) ;

% calc the body (unit) vector in the stroke plane
ahat_plane = data.AHat - strokePlaneNormals .* repmat ( dot(data.AHat, strokePlaneNormals, 2) , 1 , 3) ;
ahat_plane = ahat_plane ./ repmat(myNorm(ahat_plane),1,3) ;


% if the wing is moving "fast", then choose the chord that is more aligned
% with the wing-tip velocity

% if the wing is moving "slow" - this is the harder case. dothe following:
%
%{
find a time range around a wing flip
guess when the flip was, i.e. in the middle between two adjacent frames. if
frames are f and f+1, then assume flip was between them

before the flip chord vector should point consistently in one direction
after the flip the chord vector should point consistently in the opposite
direction. "direction" here is with respect to the tip velocity vector in
the stroke plane.

the "before" direction is the wing-tip velocity in the beginning of the
time interval

the "after" direction is the wing-tip velocity in the end of the time
interval
%}

newChord     = zeros(N, 3) ;
newAltChord  = zeros(N, 3) ;
newDiag1     = zeros(N, 1) - 1; % a negative value indicates a frame we did not get to.
newDiag2     = zeros(N, 1) - 1;

swap_flag = false(N,1) ;

vv       = dot(tip_v_plane, ahat_plane,2) ; % calc component of tip velocity along the ahat, both vectors in the stroke plane
fast_ind = ( abs(vv) > SPEED_THRESHOLD) ; % determine "fast" criterion

%fwd_stroke_flag = ( vv> 0);

% Easy case first - handle the frames with higher wing-tip velocity
for it=1:N
    if (~fast_ind(it))
        continue ;
    end
    v = tip_v_plane(it,:) ; % tip_v(it,:) ; % curret tip velocity
    v = v / norm(v) ; % normalize to get direction only
    
    c1 = mainChord(it,:) ;
    c2 = altChord(it,:) ;
    
    % if needed, invert c1 and c2 such that they point in the direction of
    % tip velocity
    
    if (dot(c1, v)<0)
        c1 = - c1 ;
    end
    if (dot(c2, v)<0)
        c2 = - c2 ;
    end
    
    %{
    % in fwd stroke, the chord should point forward
    sgn = fwd_stroke_flag(it)*1 + ~fwd_stroke_flag(it)*(-1) ; % 1 for fwd, -1 for back stroke
        
    if (sgn*dot(c1, ahat_plane(it,:))<0)
        c1 = - c1 ;
    end
    if (sgn*dot(c2, ahat_plane(it,:))<0)
        c2 = - c2 ;
    end
    %}
    % if there is only one chord vector that points up w.r.t. stroke plane,
    % choose this one
    dot1 = dot(c1, strokePlaneNormals(it,:)) ;
    dot2 = dot(c2, strokePlaneNormals(it,:)) ;
    
    if ( (dot2>0 && dot1<0) || (dot2>0 && diag2(it)/diag1(it)>1.2)) 
        swap_flag(it)     = true ;
        newChord(it,:)    = c2 ; % swap
        newAltChord(it,:) = c1 ;
        newDiag1(it)      = diag2(it) ;
        newDiag2(it)      = diag1(it) ;
    else
        newChord(it,:)    = c1 ; % do not swap
        newAltChord(it,:) = c2 ;
        newDiag1(it)      = diag1(it) ;
        newDiag2(it)      = diag2(it) ;
    end
end

% find the mean diag value of the frames processed so far
ind = (newDiag1>0) ; % values of 0 correspond to cases where wing was not identified...
meanDiag = mean(newDiag1(ind)) ;
stdDiag  = std(newDiag1(ind)) ;
% hdiag = figure;
% subplot(2,1,1) ;
% hist(newDiag1(ind),50) ;
% xlim = get(gca,'xlim') ;
% xlabel('Diag1 length') ; ylabel('counts') ; title('pass 1') ;



%% pass 2
% go over all fast frames again, and check for the ones with diag1 too
% small and diag2 around the "good size"
for it=1:N
    if (~fast_ind(it))
        continue ;
    end
    if (newDiag1(it) < 0.5*newDiag2(it)) &&  ... % if diag1 is smaller than half diag2
            abs(newDiag2(it)-meanDiag)/stdDiag < 3 % and if diag2 is with 3 sigma from mean diag
        % then swap
        swap_flag(it) = ~swap_flag(it) ;
        
        tmp = newChord(it,:) ;
        newChord(it,:) = newAltChord(it,:) ;
        newAltChord(it,:) = tmp;
        
        tmp =  newDiag1(it) ;
        newDiag1(it) = newDiag2(it) ;
        newDiag2(it) = tmp ;
        disp(it);
        clear tmp
    end
end

% check the histogram again
ind = (newDiag1>0) ; % values of 0 correspond to cases where wing was not identified...
meanDiag = mean(newDiag1(ind)) ;
stdDiag  = std(newDiag1(ind)) ;
% figure(hdiag);
% subplot(2,1,2) ;
% hist(newDiag1(ind),50,'r') ;
% xlabel('Diag1 length') ; ylabel('counts') ; title('pass 2') ;
% set(gca,'xlim',xlim) ;
% 

%% handle the frames with the slower wing-tip.
% assume the slow frames are contiguous.

% find the slow intervals
cc = bwconncomp(~fast_ind) ;

for k=1: cc.NumObjects
    t1 = cc.PixelIdxList{k}(1) ; % interval start time
    t2 = cc.PixelIdxList{k}(end) ; % interval end time
    v1 = tip_v_plane(t1,:) ; % tip velocities in start/end times
    v2 = tip_v_plane(t2,:) ;
    
    tvec_fine = (tvec(t1): (dt/100) : tvec(t2))' ;
    
    % find mean stroke plane normal vector - n
    if (t2>t1)
        n = mean( strokePlaneNormals(t1:t2,:)) ;
    else % t1==t2
        disp(['note t1==t2==' num2str(t1)]) ;
        n = strokePlaneNormals(t1,:) ;
    end
    
    n = n / norm(n) ;
    
    tip_vx_fine = fnval(sp_tipvx, tvec_fine) ;
    tip_vy_fine = fnval(sp_tipvy, tvec_fine) ;
    tip_vz_fine = fnval(sp_tipvz, tvec_fine) ;
    
    % find the wing-tip velocity component in the mean stroke-plane (fine time intervals).
    tip_v_fine       = [ tip_vx_fine, tip_vy_fine, tip_vz_fine] ;
    nmat             = repmat (n, length(tvec_fine), 1) ;
    try
    tip_v_fine_plane = tip_v_fine - nmat .* repmat ( dot(tip_v_fine, nmat, 2) , 1 , 3) ;
    catch
        keyboard ;
    end
    speed_fine_plane = myNorm(tip_v_fine_plane) ;
    
    [~, indmin] = min(speed_fine_plane) ;
    
    % the time at which the minimum speed occurs. assume this is the flip time (!)
    tflip = tvec_fine(indmin) ;
    
    
    % go over the times in the current interval and decide if to swap the
    % chord or not.
    for it = t1:t2
        if (tvec(it)<=tflip)
            v_use = v1 ;
        else
            v_use = v2 ;
        end
        % the chord should align more with v_use
        
        c1 = mainChord(it,:) ;
        c2 = altChord(it,:) ;
        
        dot1 = dot(c1, v_use) ;
        dot2 = dot(c2, v_use) ;
        
        continue_process = true ;
        
        if ( abs(tvec(it)-tflip)*data.params.fps <= 0.25) % if very close to the flipping point <0.25 frames
            % take the longer diag
            if (diag2(it)/diag1(it)>1.2)
                swap_flag(it) = true ;
            end
            continue_process = false ;
        end
        
        if (continue_process)            
            % if diag criterion applies, swap anyway, ignoring velocity criterion
            if (diag1(it) < 0.66*diag2(it)) &&  ... % if diag1 is smaller than half diag2
                    abs(diag2(it)-meanDiag)/stdDiag < 3   % and if diag2 is with 3 sigma from mean diag
                swap_flag(it) = true ;
                continue_process = false ;
            end
        end
        
        if (continue_process)            
            % if diag criterion applies, swap anyway, ignoring velocity criterion
            if (diag2(it) < 0.66*diag1(it)) &&  ... % if diag1 is smaller than half diag2
                    abs(diag1(it)-meanDiag)/stdDiag < 3   % and if diag2 is with 3 sigma from mean diag
                swap_flag(it) = false ;
                continue_process = false ;
            end
        end
        
        if (continue_process)    
            if (dot2>0 && dot1<0)
                swap_flag(it) = true;
            end
        end
        
        %{
        % older code
        diagRatio = diag2(it) / diag1(it) ;
        if (diagRatio>1.5) % if the diag ratio is too large, swap anyway, ignoring speed information
            swap_flag(it) = true ;
        elseif (dot2>dot1 && diagRatio>0.8) % velocity needs to be along the correct direction, but diag2 can't be too small
            swap_flag(it) = true ;
        end
        
        if (swap_flag(it))
            newChord(it,:) = c2 ; % swap
        else
            newChord(it,:) = c1 ; % do not swap
        end
        %}
        
        if (swap_flag(it))
            swap_flag(it)     = true ;
            newChord(it,:)    = c2 ; % swap
            newAltChord(it,:) = c1 ;
            newDiag1(it)      = diag2(it) ;
            newDiag2(it)      = diag1(it) ;
        else
            newChord(it,:)    = c1 ; % do not swap
            newAltChord(it,:) = c2 ;
            newDiag1(it)      = diag1(it) ;
            newDiag2(it)      = diag2(it) ;
        end
    end
end

% last pass, go over the "fast" frames again, and check for a continuity
% criterion for the "pitch" angle

dot1 = dot(newChord, strokePlaneNormals, 2) ;
dot2 = dot(newAltChord, strokePlaneNormals, 2) ;

ang1 = acos(dot1) * 180 / pi;
ang2 = acos(dot2) * 180 / pi ;
% 
% figure ; hold on ;
% plot(dot1,'bo-') ;
% plot(dot2,'rs-') ;
% legend({'dot1','dot2'}) ; grid on ; box on ;
%keyboard ;

for it=2:N-1
    if (~fast_ind(it))
        continue ;
    end
    
    % CODE UNDER CONSTRUCTION
    % STOPPED HERE! 
    
    
end




% PLOTS
if (plotFlag)
    figure;
    subplot(2,3,1) ; hold on ;
    xx = 1:N ; % tvec_ms ;
    
    plot(xx, tip(:,1),'bo') ;
    plot(xx, tipx_smooth,'k-','linewidth',2) ;
    plot(xx, tipx_smooth+ESTERR,'k--') ;
    plot(xx, tipx_smooth-ESTERR,'k--') ;
    title('wing tip x') ;
    
    subplot(2,3,2) ; hold on ;
    plot(xx, tip(:,2),'bo') ;
    plot(xx, tipy_smooth,'k-','linewidth',2) ;
    plot(xx, tipy_smooth+ESTERR,'k--') ;
    plot(xx, tipy_smooth-ESTERR,'k--') ;
    title('wing tip y') ;
    
    subplot(2,3,3) ; hold on ;
    plot(xx, tip(:,3),'bo') ;
    plot(xx, tipz_smooth,'k-','linewidth',2) ;
    plot(xx, tipz_smooth+ESTERR,'k--') ;
    plot(xx, tipz_smooth-ESTERR,'k--') ;
    title('wing tip z') ;
    
    subplot(2,3,4) ;
    plot(xx, tip_vx,'r-','linewidth',2) ;
    title('wing tip v_x') ;
    
    subplot(2,3,5) ;
    plot(xx, tip_vy,'r-','linewidth',2) ;
    title('wing tip v_y') ;
    
    subplot(2,3,6) ;
    plot(xx, tip_vz,'r-','linewidth',2) ;
    title('wing tip v_z') ;
    
    
    figure ;
    plot(xx, speed,'m-^','linewidth',2) ;
    title('wing-tip speed w.r.t. body CM') ;
    ylabel('speed [m/s]') ; xlabel('time [ms]') ;
    grid on ;
    
    % generate a plot for every wing-stroke (groups of 2*DELTA frames)
    startInd = 1:(2*DELTA):N ;
    endInd   = [startInd(2:end), N] ;
    for f=1:length(startInd) ;
        t1 = startInd(f) ;
        t2 = endInd(f) ;
        ind = t1:t2 ;
        figure; hold on ;
        plot3(tip(ind,1), tip(ind,2), tip(ind,3),'ko') ;
        
        plot3(othertip(ind,1), othertip(ind,2), othertip(ind,3) ,'ko') ;
        
        plot3(tipx_smooth(ind), tipy_smooth(ind), tipz_smooth(ind),'b.-') ;
        plot3(tip(t1,1), tip(t1,2), tip(t1,3),'ks','markerfacecolor','r') ;
        
        a = data.AHat(t1,:) * 40 * data.params.voxelSize ;
        
        plot3( [1 -1]*a(1) , [1 -1]*a(2) , [1 -1]*a(3) ,'-','linewidth',3,'color',[0 0.7 0]) ;
        plot3(0,0,0,'o','color',[0 0.7 0],'linewidth',3) ;
        C = 0.7e-3 ;
        for it=t1:t2
            col = [0, 0, 0]*swap_flag(it) + [1 0 0]*(~swap_flag(it)) ;
            if (~fast_ind(it) && ~swap_flag(it)) 
                col = [1 0.5 0] ;
            end
            
            plot3( tip(it,1)+[0, newChord(it,1)]*C, ...
                tip(it,2)+[0, newChord(it,2)]*C, ...
                tip(it,3)+[0, newChord(it,3)]*C, '-','linewidth',2,'color',col) ;
            text( tip(it,1), tip(it,2), tip(it,3),['  ' num2str(it)]) ;
        end
        
        axis equal ; box on ; grid on ; view(-51, 22) ; rotate3d on ;
        axis tight;
        title(['Frame indices ' num2str(t1) ' - ' num2str(t2)]) ;
    end
   
end

%keyboard ;
end


%{
% smooth using spaps that allows to specify tolerance "tol"
% note that if the estimated error of the data is +/-(err),
% and the time interval between samples is dt then
% weight = dt
% roughly tol=(2*err)^2 * dt
dt = tplotSec(2)-tplotSec(1) ;

esterr = 1 ; % measurement is (+/-)esterr about the "real" data [in degrees]
N      = length(tplotSec) ;

tol = N * esterr^2 * dt

splinePhi   = spaps(tplotSec, phiBodyDeg, tol) ;
%}

%%
% this function was taken from the MathWorks website
% http://www.mathworks.com/matlabcentral/fileexchange/43305-plane-fit/content/affine_fit.m
function [n,V,p] = affine_fit(X)
%Computes the plane that fits best (lest square of the normal distance
%to the plane) a set of sample points.
%INPUTS:
%
%X: a N by 3 matrix where each line is a sample point
%
%OUTPUTS:
%
%n : a unit (column) vector normal to the plane
%V : a 3 by 2 matrix. The columns of V form an orthonormal basis of the
%plane
%p : a point belonging to the plane
%
%NB: this code actually works in any dimension (2,3,4,...)
%Author: Adrien Leygue
%Date: August 30 2013

%the mean of the samples belongs to the plane
p = mean(X,1);

%The samples are reduced:
R = bsxfun(@minus,X,p);
%Computation of the principal directions if the samples cloud
[V,D] = eig(R'*R);
%Extract the output from the eigenvectors
n = V(:,1);
V = V(:,2:end);
end