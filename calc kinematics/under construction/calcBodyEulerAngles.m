function [eulerAngles smoothedRollAngle fdLambda ci tvec rhoTimes] = calcBodyEulerAngles (data, rhoTimes, plotFlag)

plotFlag2 = false ;

rollVectors = data.rollVectors ; % copy

% make sure the first and last frames are included. if not - copy 

if (~isfield(data,'startAnalysisTimeMS'))
    if (rhoTimes(1)~=1)
        rollVectors(1,:) = rollVectors(rhoTimes(1),:) ;
        rhoTimes = [1 rhoTimes] ;
    end
    if (rhoTimes(end)~=data.Nimages)
        rollVectors(end,:) = rollVectors(rhoTimes(end),:) ;
        rhoTimes = [rhoTimes data.Nimages] ;
    end
end
if (~exist('plotFlag','var'))
    plotFlag = true;
end

if (isfield(data,'startAnalysisTimeMS'))
    %startTime = data.startAnalysisTimeMS * data.params.fps / 1000 ;
    %endTime = data.endAnalysisTimeMS * data.params.fps / 1000 ;
    Np = endTime - startTime + 1 ;
else
    startTime = data.params.startTrackingTime ;
    endTime   = data.params.endTrackingTime ;
    Np = data.Nimages ;
end

%N = data.Nimages ;
%tvec = (startTime:endTime) / data.params.fps * 1000 ; % in ms
%tvec = rhoTimes 
%tvec = tvec' ;
startFrame = startTime - data.params.startTrackingTime + 1 ;
endFrame   = startFrame + Np - 1 ;
frameIndices = startFrame:endFrame ;

Nt  = length(rhoTimes) ;
eulerAngles = zeros(Nt,3) ;  % yaw, pitch, roll

if (plotFlag2)
    h = figure;
    orig = [0 0 0 ]; 
else
    h = -1 ;
end

C = 180 / pi ;

for j=1:Nt
    t = rhoTimes(j) ;
    xb = data.AHat(t, :) ;
    yb = rollVectors(t,:) ; %  goes from R-->L wing hinges
    zb = cross(xb,yb) ;
    
    %R = [xb' yb' zb'] ;    
    %roll = atan2(R(3,1), R(3,2)) * 180 / pi 
    %pitch = acos(R(3,3)) * 180 / pi 
    %yaw   = -atan2(R(1,3), R(2,3)) * 180 / pi 
    
    yaw   = atan2(xb(2), xb(1)) * C ;
    pitch = asin(xb(3)) * C ;
    
    xb_proj = [xb(1) xb(2) 0 ] ;
    xb_proj = xb_proj / norm(xb_proj) ;
    v = cross([0 0 1], xb_proj) ;
    roll = acos( dot(v, yb)) * C * sign(yb(3));
    
    eulerAngles(j,:) = [yaw pitch roll] ;
    
    if (plotFlag2)
       figure(h) ; clf ; hold on ;
       myplot(orig, xb,'go-');
       myplot(orig, yb,'ro-') ;
       myplot(orig, zb,'bs-') ;
       myplot(orig, xb_proj,'go--') ;
       myplot(orig, v, 'k--') ;
       view(3) ; grid on ; box on
       xlabel('x') ; ylabel('y') ; zlabel('z');    
       title(t) ;
       pause(0.4) ;
    end
end

if (nargout>1)
    fps = data.params.fps ;
    %tvec = ( rhoTimes - 1 + data.params.startTrackingTime ) / fps * 1000 ; % in ms
    tvec = ( rhoTimes - 1 + data.params.startTrackingTime ) / fps  ; % in SECONDS
    tvec = tvec' ; 
    
    smoothRhoFlag = ~isfield(data,'rho_fd_lambda') ;
    
    if (~smoothRhoFlag)
        smoothRhoFlag = isempty(data.rho_fd_lambda) ;
    end
    %smoothRhoFlag=true ;
    if (smoothRhoFlag)
        disp('start smoothing rho') ;
        disp(['calcBodyEulerAngles: smoothing rho(t) using [order, rmserr, LFD] = [' ...
            num2str(data.rho_norder) ', ' num2str(data.rho_rmserr) ', ' ...
            num2str(data.rho_LFD_order) ']']) ;
        norder = data.rho_norder ; % 5;   % ; 6
        RMSERR = data.rho_rmserr ; % 0.4; % 1.5
        LFD_ORDER = data.rho_LFD_order ; % 2 ; % 2 ;
        %plotFlag = true ;
        lambdaHi = 100 ; % 100000 ;
        lambdaLo = 1e-20 ; 
        lambdaHiLo = [lambdaHi lambdaLo] ;
        %nbasisoverride = round(length(tvec)/2) ;
        
        [fdLambda lambdaMid rmsMid ci] = ...
            mySmoothWrapperMSE(tvec, eulerAngles(:,3), norder,RMSERR, ...
            LFD_ORDER, [],plotFlag, lambdaHiLo);%#ok<ASGLU> %, nbasisoverride) ; %#ok<ASGLU>
        
        %[fdLambda lambdaMid rmsMid ci] = ...
        %    mySmoothWrapperMSE_mk2(tvec, eulerAngles(:,3), norder,RMSERR, ...
        %   LFD_ORDER, [], plotFlag,lambdaHiLo);
    else
        fdLambda = data.rho_fd_lambda ;
        lambdaMid = [] ;
        rmdMid = [] ;
        ci = [] ;
    end
    
    allT = (0:data.Nimages-1)' + data.params.startTrackingTime ;
    allT = allT / fps  ; % in SEC
    
    subT = (startTime:endTime)' / fps   ; % in SECONDS
    subRho = eval_fd(subT, fdLambda) ;
    
    smoothedRollAngle = zeros(data.Nimages,1) ;
    smoothedRollAngle(frameIndices) = subRho ;
        
    if (plotFlag)
        figure ; hold on ;
        plot(tvec*1000, eulerAngles(:,3),'ko','markerfacecolor',[0 0.8 0]) ;
        plot(subT*1000, subRho,'k-');
        set(gca,'fontsize',14) ;
        xlabel('t [ms]') ;
        ylabel('\rho [deg]') ;
        grid on ; box on ;
        set(gca,'xlim',[subT(1) subT(end)]*1000);
    end
else
    smoothedRollAngle = [] ;
end

end

function myplot(v1,v2,colstr)
v3 = v1+v2 ;
plot3( [v1(1) v3(1)], [v1(2) v3(2)], [v1(3) v3(3)],colstr,'linewidth',2) ;
end