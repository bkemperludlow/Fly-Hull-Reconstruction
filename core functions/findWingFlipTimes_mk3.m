function [fwdFlipTimes, backFlipTimes, fwdFlipInd, backFlipInd, fwdFlipPhi, ...
    backFlipPhi, badIndices] = findWingFlipTimes_mk3 (t, phi, plotFlag)
% fwd flips is when phi is minimum (towards 0)
% back flips is when phi is maximum (towards 180)

if (~exist('plotFlag','var'))
    plotFlag = false ;
end

W = 9 ;% windows is it-9 : it+9
N = length(t) ;

maxIndices = zeros(N,1) ;
minIndices = zeros(N,1) ;
badIndices = zeros(N,1) ;
maxCounter = 0 ;
minCounter = 0 ;
badCounter = 0 ;
testvec = NaN(N,1) ;

% find neighborhood of local extrema
for it=2:(N-1)
    i1 = max([1, it-W]) ;
    i2 = min([N, it+W]) ;
    
    ff = phi(i1:i2) ;
    
    % filter out single-point outliars
    stdvec = zeros(length(ff),1) ;
    for k=1:length(ff)
        % exclude the k'th point
        ff2 = ff ;
        ff2(k) = [] ;
        % re-calc std
        d = diff(ff2) ;
        stdvec(k) = std(d) / mean(abs(d)) ;
        if (k == it-i1+1 )
            testvec(it) = stdvec(k) ;
        end
    end
    
end

for it=2:(N-1)
    if (testvec(it-1)-testvec(it)>0.2) && (testvec(it+1)-testvec(it)>0.2)
        badCounter = badCounter + 1 ;
        badIndices(badCounter) = it;
    end
end

badIndices      = badIndices(1:badCounter) ;
origPhi         = phi ;
phi(badIndices) = NaN ;

    %e     = 
    % find if there are differences that are too large in terms of std/mean

for it=2:(N-1)
    i1 = max([1, it-W]) ;
    i2 = min([N, it+W]) ;
    
    ff = phi(i1:i2) ;    
  
    mx = max(ff) ; 
    % use indmax and indmin in case there are two max elements
    % if there is more than one max, use only the first one.
    indmax = find(ff==mx) ;
    if (numel(indmax)>1)
        ff(indmax(2:end))=NaN ;
    end
    mn = min(ff) ; 
    indmin = find(ff==mn) ;
    if (numel(indmin)>1)
        ff(indmin(2:end))=NaN ;
    end
        
    
    oneNotNaN = xor(isnan(phi(it)), isnan(ff(it-i1+1))) ; 
    bothNotNaN = and (~isnan(phi(it)), ~isnan(ff(it-i1+1))) ; 
    if (oneNotNaN || ( bothNotNaN && (phi(it)~=ff(it-i1+1)))) 
        disp('problem')
        keyboard ;
    end
    
    if (ff(it-i1+1) == mx)
        maxCounter = maxCounter + 1 ;
        maxIndices(maxCounter) = it ;
    end
    
    if (ff(it-i1+1) == mn)
        minCounter = minCounter + 1 ;
        minIndices(minCounter) = it ;
    end
    
end

% trim
maxIndices = maxIndices(1:maxCounter) ;
minIndices = minIndices(1:minCounter) ;


fwdFlipTimes  = t(minIndices) ;
backFlipTimes = t(maxIndices) ;
fwdFlipInd    = minIndices ;
backFlipInd   = maxIndices ;
fwdFlipPhi    = phi(minIndices) ;
backFlipPhi   = phi(maxIndices) ;


if plotFlag
    figure ; hold on ;
    plot(origPhi, 'o-','color',[0.1 0.7 0.1]) ;
    plot(maxIndices, origPhi(maxIndices),'k^','markerfacecolor','r') ;
    plot(minIndices, origPhi(minIndices),'kv','markerfacecolor','b') ;
    plot(badIndices, origPhi(badIndices),'kx','linewidth',3,'markersize',17) ;
    set(gca,'fontsize',14) ;
    grid on ; box on ;
    set(gca,'xlim',[1 N]) ; xlabel('Time [frames]') ; ylabel('Angle [deg]') ;
    legend({'\phi','Fwd flips','Back flips','Rejected'},'fontsize',10,'location','northwest');
end