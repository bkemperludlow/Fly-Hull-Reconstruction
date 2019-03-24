function [ fwdFlipIndR, backFlipIndR, fwdFlipTimesR,  backFlipTimesR,  ...
    fwdFlipIndL,  backFlipIndL,  fwdFlipTimesL, backFlipTimesL ] = ...
    findWingFlipTimes_mk2 (sp_phiR, sp_phiL, t)

%disp('this function has some issues with double peaks in the smoothed phi. See use of Np.') ;

% input arguments: the fd spline structures of phi right and phi left

%fps = 8000 ;
%t = (0:data.Nimages-1) + data.params.startTrackingTime ;
%t = t / fps * 1000 ; % in ms
trange = [t(1) t(end)] ;
dt = (t(2) - t(1)) / 10  ;

Np = 10 * 8 / 2 ; % 0.5 ms in 8000fps, 

% define a denser time grid for plotting smoothed data
t2 = trange(1):dt:trange(end) ;

Nalloc = 100 ; % number of flips we allocate. we trim later on.

% ---------------------------------
% FIND THE FLIP TIMES FOR EACH WING
% ---------------------------------

% note there was a confusion between forward and backward flips. instead of
% changing variable names throughout the algorithm, I only swapped the
% final resuls.

phiR0 = fnval(sp_phiR, t2); %-eval_fd(t2, fdPhiR) ;
phiR1 = fnval(fnder(sp_phiR,1), t2) ; % eval_fd(t2, fdPhiR,1) ;

phiL0 = fnval(sp_phiL, t2); % eval_fd(t2, fdPhiL) ;
phiL1 = fnval(fnder(sp_phiL,1), t2) ; %eval_fd(t2, fdPhiL,1) ;

Nt = length(t2) ;
backFlipTimesL = zeros(1,Nalloc) ;
backFlipTimesR = zeros(1,Nalloc) ;
backFlipIndR   = zeros(1,Nalloc) ;
backFlipIndL   = zeros(1,Nalloc) ;
backRcounter   = 0 ;
backLcounter   = 0 ;

fwdFlipTimesL = zeros(1,Nalloc) ;
fwdFlipTimesR = zeros(1,Nalloc) ;
fwdFlipIndR   = zeros(1,Nalloc) ;
fwdFlipIndL   = zeros(1,Nalloc) ;
fwdRcounter   = 0 ;
fwdLcounter   = 0 ;

for it=1:Nt-1
    y1=phiR1(it);
    y2=phiR1(it+1);
    ang=phiR0(it+1) ;
    ang0=phiR0(it) ;
    % right back flip
    if (y2*y1<0 && y2>y1 && ang<inf)
        %if (it>Np && it<Nt-Np && mean(phiR0(it-Np:it-1)-ang0)*mean(phiR0(it+1:it+Np)-ang0)>0)
            dy = y2-y1 ;
            slope = dy / dt ;
            backRcounter = backRcounter + 1 ;
            backFlipTimesR(backRcounter) = t2(it) + y1/slope ;
            backFlipIndR(backRcounter) = it ;
        %end
    end
    % right fwd flip
    if (y2*y1<0 && y2<y1 && ang<inf)
        %if (it>Np && it<Nt-Np && mean(phiR0(it-Np:it-1)-ang0)*mean(phiR0(it+1:it+Np)-ang0)>0)
            dy = y2-y1 ;
            slope = - dy / dt ;
            fwdRcounter = fwdRcounter + 1 ;
            fwdFlipTimesR(fwdRcounter) = t2(it) + y1/slope ;
            fwdFlipIndR(fwdRcounter) = it ;
        %end
    end
    
    y1=phiL1(it);
    y2=phiL1(it+1);
    ang=phiL0(it+1) ;
    % left back flip
    if (y2*y1<0 && y2>y1 && ang<inf)
        dy = y2-y1 ;
        slope = dy / dt ;
        backLcounter = backLcounter + 1 ;
        backFlipTimesL(backLcounter) = t2(it) + y1/slope ;
        backFlipIndL(backLcounter) = it ;
    end
    % fwd left flip
    if (y2*y1<0 && y2<y1 && ang<inf)
        dy = y2-y1 ;
        slope = - dy / dt ;
        fwdLcounter = fwdLcounter + 1 ;
        fwdFlipTimesL(fwdLcounter) = t2(it) + y1/slope ;
        fwdFlipIndL(fwdLcounter) = it ;
    end
end

% trim
fwdFlipIndR = fwdFlipIndR(1:fwdRcounter) ;
fwdFlipIndL = fwdFlipIndL(1:fwdLcounter) ;

fwdFlipTimesR = fwdFlipTimesR(1:fwdRcounter) ;
fwdFlipTimesL = fwdFlipTimesL(1:fwdLcounter) ;

backFlipIndR = backFlipIndR(1:backRcounter) ;
backFlipIndL = backFlipIndL(1:backLcounter) ;

backFlipTimesR = backFlipTimesR(1:backRcounter) ;
backFlipTimesL = backFlipTimesL(1:backLcounter) ;

% swap forward and backward wing-flip times
[fwdFlipIndR,    backFlipIndR   ] = mySwap(fwdFlipIndR,    backFlipIndR ) ;
[fwdFlipTimesR,  backFlipTimesR ] = mySwap(fwdFlipTimesR,  backFlipTimesR) ;
[fwdFlipIndL,   backFlipIndL  ]   = mySwap(fwdFlipIndL,    backFlipIndL) ;
[fwdFlipTimesL, backFlipTimesL]   = mySwap(fwdFlipTimesL,  backFlipTimesL) ;

% take out wing flips that are too close apart
%flags = diff(fwdFlipTimesR) < DELTA

end

